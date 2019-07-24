#include "Learner.h"

#include <cmath>
#include <cfloat>
#include <algorithm>
#include <numeric>
#include <vector>
#include <fstream>

#include "mkl.h"

#include "Recognize.h"

Learner::Learner()
    : common_support(0.0), plogp(0.0), kl(0.0)
{
}

Learner::~Learner()
{
}

void Learner::Renormalize()
{
    const MKL_INT n = GetNumberOfParameters();
    const MKL_INT k = GetNumberOfConstraints();

    aux.resize(n + k);
    double* y = aux.data();
    double* g = aux.data() + n;

    // y <- exp(x)
    vdExp(n, _x.data(), y);
    
    // g <- C^t.y
    C.dot(y, g, SPARSE_OPERATION_TRANSPOSE);
    
    // g <- log(C^t.exp(x))
    vdLn(k, g, g);

    // x -= C.g
    C.dot(g, _x.data(), SPARSE_OPERATION_NON_TRANSPOSE, -1.0, 1.0);
}

void Learner::RewriteWeights(Fsa & fsa) const
{
    for (auto& t : fsa.GetTransitionMtx())
    {
        auto& transitions = t.second.transitions;
        auto& emissions = t.second.emissions;

        for (auto& edge : emissions)
            edge.logprob = (edge.index >= 0 ? GetWeight(edge.index) : 0.0);

        for (auto& edge : transitions)
            edge.logprob = (edge.index >= 0 ? GetWeight(edge.index) : 0.0);
    }
}

void Learner::PrintJ(FILE* f)const 
{
    //vdExp(_x.size(), _x.data(), _x.data());
    PrintCsrMtx(f, _x, Crow, Ccol);
    //vdLn(_x.size(), _x.data(), _x.data());
}

void Learner::PrintC(FILE* f)const
{
    PrintCsrMtx(f, ones, Crow, Ccol);
}

void Learner::PrintP(FILE* f) const
{
    PrintCsrMtx(f, Pdata, Prow, Pcol);
}

void Learner::PrintM(FILE* f)const
{
    PrintCsrMtx(f, ones, Mrow, Mcol);
}

bool Learner::SaveMatrices(const std::string& filename) const
{
    std::ofstream ofs(filename + ".C");
    ofs.precision(DBL_DIG);
    if (ofs)
    {
        WriteCsrMtx(ofs, ones, Crow, Ccol);
        ofs.close();
    }
    else
        return false;
    ofs.open(filename + ".M");
    if (ofs)
    {
        WriteCsrMtx(ofs, ones, Mrow, Mcol);
        ofs.close();
    }
    else
        return false;
    ofs.open(filename + ".P");
    if (ofs)
    {
        WriteCsrMtx(ofs, Pdata, Prow, Pcol);
        ofs.close();
    }
    else
        return false;
    ofs.open(filename + ".prob");
    if (ofs)
    {
        for (const auto& prob : p)
            ofs << prob << std::endl;
        ofs.close();
    }
    else
        return false;
    ofs.open(filename + ".aux");
    if (ofs)
    {
        ofs << common_support << std::endl;
        ofs << plogp << std::endl;
        ofs << model_volume << std::endl;
        ofs << aux_hessian << std::endl;
        ofs << auxiliary_parameters << std::endl;

        ofs.close();
    }
    else
        return false;
    return true;
}

bool Learner::LoadMatrices(const std::string& filename)
{
    std::ifstream ifs(filename + ".C");
    if (ifs)
    {
        ReadCsrMtx(ifs, ones, Crow, Ccol);
        ifs.close();
    }
    else
        throw LearnerError("Unable to open \"", filename + ".C", "\"!");

    ifs.open(filename + ".M");
    if (ifs)
    {
        ReadCsrMtx(ifs, ones, Mrow, Mcol);
        ifs.close();
    }
    else
        throw LearnerError("Unable to open \"", filename + ".M", "\"!");

    ifs.open(filename + ".P");
    if (ifs)
    {
        ReadCsrMtx(ifs, Pdata, Prow, Pcol);
        ifs.close();
    }
    else
        throw LearnerError("Unable to open \"", filename + ".P", "\"!");

    ifs.open(filename + ".prob");
    do
    {
        p.emplace_back();
        ifs >> p.back();
    } while (ifs);
    p.pop_back();
    ifs.close();

    ifs.open(filename + ".aux");
    if (ifs)
    {
        ifs >> common_support;
        ifs >> plogp;
        ifs >> model_volume;
        ifs >> aux_hessian;
        ifs >> auxiliary_parameters;
        if (!ifs)
            throw LearnerError("Invalid data in \"", filename + ".aux", "\"!");
        ifs.close();
    }
    else
        throw LearnerError("Unable to open \"", filename + ".aux", "\"!");

    if (p.size() + 1 != Mrow.size())
        throw LearnerError("Size mismatch: ", p.size() + 1, " != ", Mrow.size());

    if (*std::max_element(Pcol.begin(), Pcol.end()) >= GetNumberOfParameters())
        throw LearnerError("Size mismatch: more path indexes than parameters in the automaton!");

    if (Mcol.back() + 1 != Prow.size() - 1)
        throw LearnerError("M cols != P rows");

    _x.assign(GetNumberOfParameters(), 0.0);

    return true;
}

std::vector<double> Learner::GetOptimizationInfo()
{
    return std::vector<double>();
}

std::string Learner::GetOptimizationHeader()const
{
    return std::string();
}

std::vector<double> Learner::GetOptimizationResult(bool)
{
    return std::vector<double>();
}

bool Learner::HaltCondition(double)
{
    return false;
}

template<class Indexed>
void Collect(MKL_INT& k, double& model_volume,
    std::vector<double>& _x,
    std::vector<MKL_INT>& Crow, std::vector<MKL_INT>& Ccol,
    const Indexed& emissions)
{
    if (emissions.size() > 1)
    {
        auto start_index = Ccol.size();
        for (const auto& next : emissions)
        {
            if (next.index != (MKL_INT)Ccol.size())
                throw LearnerError("Indexing error: ", next.index, " != ", Ccol.size());

            Crow.push_back(Ccol.size());
            Ccol.emplace_back(k);
            _x[next.index] = next.logprob;
        }
        ++k; // number of constraints
        const auto dim = Ccol.size() - start_index; // number of parameters in this constraint
        model_volume += LogSimplexVolume(dim);
    }
}

void Learner::BuildConstraints(const Fsa& fsa)
{
    //TODO orthogonal basis for constraints!
    model_volume = 0;
    Crow.clear(); Ccol.clear();

    _x.assign(fsa.GetNumberOfParameters(), 0.0);
    Crow.reserve(fsa.GetNumberOfParameters()); Ccol.reserve(fsa.GetNumberOfParameters());
    MKL_INT k = 0;

    for (const auto& t : fsa.GetTransitionMtx())
    {
        const auto& transitions = t.second.transitions;
        const auto& emissions = t.second.emissions;

        // collect equivocal emissions, then transitions
        Collect(k, model_volume, _x, Crow, Ccol, emissions);
        Collect(k, model_volume, _x, Crow, Ccol, transitions);
    }
    Crow.push_back(Ccol.size());
}

typedef std::vector<std::pair<MKL_INT, double>> Path;

void Learner::BuildFrom(const Fsa& fsa, const Corpus& corpus, bool bfs)
{
    BuildConstraints(fsa);
    BuildPaths(fsa, corpus, bfs);
    Trim();
}

void Learner::BuildPaths(const Fsa& fsa, const Corpus& corpus, bool bfs)
{
    trimmed_weights.assign(GetNumberOfParameters(), -2); // per default every index is unused
    auxiliary_parameters = 0;
    const auto& startState = fsa.GetTransitionMtx().at(fsa.GetStartState());
    const auto endstate = fsa.GetEndState();
    bool has_path = false;
    Corpus::const_iterator wordIt; // currently processed word

    Recognizer<Path>::Accumulator accumulator(
        [&](Path& history, const Fsa::Transitions::value_type& transition, const Fsa::Emissions::value_type& emission) -> Path&
    {
        if (transition.index >= 0)
            SortedInsert(history, transition.index) += 1;
        if (!StrEq::streq(transition.next->first, endstate) && emission.index >= 0)
            SortedInsert(history, emission.index) += 1;
        return history;
    });

    Recognizer<Path>::ResultHandler  resultHandler([&](const Path& path)
    {
        if (!has_path) // first recognition of this word
        {
            Mrow.push_back(Mcol.size());
            common_support += wordIt->second;
            p.emplace_back(wordIt->second);
        }
        has_path = true;

        Prow.push_back(Pcol.size());
        for (const auto& variable : path)
        {
            Pdata.emplace_back(variable.second);
            Pcol.emplace_back(variable.first);
            trimmed_weights[variable.first] = 0; // variable is marked used
        }

        Mcol.push_back(Mcol.size());        
    });
    Recognizer<Path> recognizer(fsa.GetEndState(), accumulator, resultHandler);
    auto recognize = bfs ? &Recognizer<Path>::RecognizeBFS : &Recognizer<Path>::RecognizeDFS;

    Pdata.clear(); Pcol.clear(); Prow.clear();
    Mcol.clear(); Mrow.clear();
 
    common_support = 0.0;
    p.clear();
    p.reserve(corpus.size()); // just an estimate

    size_t n_word = 0;
    const char* current_word = "";
    ProgressIndicator(n_word, &n_word, 100.0 / corpus.size(),
        "\rRecognize: %4.3g%% \"%-20.20s\"      ",
        [&]()
    {
        aux_hessian = 0;
        for (wordIt = corpus.begin(); wordIt != corpus.end(); ++wordIt)
        {
            ++n_word;
            has_path = false;
            current_word = wordIt->first.c_str();
            (recognizer.*recognize)(current_word, startState, Path());
            if (!has_path)
            {   // word was not recognized
                ++auxiliary_parameters;
                aux_hessian -= std::log(wordIt->second);
            }
        }
    }, current_word);

    Mrow.push_back(Mcol.size());
    Prow.push_back(Pcol.size());
}

void Learner::Trim()
{
    // new -1 may appear if a parameter in a constraint appears to be alone
    MKL_INT c = -1;
    MKL_INT nnz_in_c = -1;
    for (MKL_INT i = 0; i < GetNumberOfParameters(); ++i)
    {
        const MKL_INT this_c = Ccol[i];
        if (this_c != c)
        {   // finish this constraint
            c = this_c;
            if (nnz_in_c >= 0)
                trimmed_weights[nnz_in_c] = -1;
            nnz_in_c = -1;
        }
        if (trimmed_weights[i] >= 0)
        {
            if (nnz_in_c == -1)
                nnz_in_c = i;
            else
                nnz_in_c = -2;
        }
    }
    if (nnz_in_c >= 0)
        trimmed_weights[nnz_in_c] = -1;
    {
    decltype(Ccol) Ccol_new; Ccol_new.reserve(GetNumberOfParameters());

    // renumbering parameters, reconstruct C matrix
    MKL_INT good_indices = 0;
    MKL_INT good_constraints = -1;
    c = -1;
    for (MKL_INT i = 0; i < GetNumberOfParameters(); ++i)
    {
        if (trimmed_weights[i] == 0)
        {
            trimmed_weights[i] = good_indices++;
            if (trimmed_weights[i] < i)
                _x[trimmed_weights[i]] = _x[i];
            if (c < Ccol[i])
            {
                ++good_constraints;
                c = Ccol[i];
            }
            Ccol_new.push_back(good_constraints);
        }
    }
    std::swap(Ccol, Ccol_new);
    Crow.resize(good_indices + 1);
    } // destruct swapped Ccol

    _x.resize(GetNumberOfParameters());
    
    // reconstruct P matrix
    decltype(Pcol) Pcol_new; Pcol.reserve(Pcol.size());
    decltype(Pdata) Pdata_new; Pdata.reserve(Pdata.size());
    
    MKL_INT thisrowstart = Prow[0];
    for (MKL_INT pathindex = 0; pathindex < GetNumberOfPaths(); ++pathindex)
    {
        MKL_INT thisrowend = Prow[pathindex + 1];
        for (MKL_INT j = thisrowstart; j < thisrowend; ++j)
        {
            if (trimmed_weights[Pcol[j]] >= 0)
            {
                Pcol_new.push_back(trimmed_weights[Pcol[j]]);
                Pdata_new.push_back(Pdata[j]);
            }
        }
        thisrowstart = Prow[pathindex + 1];
        Prow[pathindex + 1] = Pcol_new.size();
    }

    std::swap(Pcol_new, Pcol);
    std::swap(Pdata_new, Pdata);    
}

double Learner::GetWeight(MKL_INT i) const
{
    static const double minf = atof("-inf");
    switch (trimmed_weights[i])
    {
    case -2: return minf;
    case -1: return 0.0;
    default: return _x[trimmed_weights[i]];
    }
}

void Learner::LambdaUpdate(double* lstep, double* l, double eta, bool exponential) const
{
    if (!exponential)
    {
        cblas_daxpy(GetNumberOfConstraints(), -eta, lstep, 1, l, 1);
    }
    else
    {
        // exponential update
        /*
        lambda *= e^{-eta*lstep/lambda}
        */

        // laux /= lambda
        vdDiv(GetNumberOfConstraints(), lstep, l, lstep);

        cblas_dscal(GetNumberOfConstraints(), -eta, lstep, 1);

        // e^laux
        vdExp(GetNumberOfConstraints(), lstep, lstep);

        // lambda *= laux
        vdMul(GetNumberOfConstraints(), l, lstep, l);
    }
}

void Learner::FinalizeCallback() {}

void Learner::Finalize()
{
    aux.resize(GetNumberOfStrings());

    ones.assign(std::max(GetNumberOfPaths(), GetNumberOfParameters()), 1.0);

    vdLn(GetNumberOfStrings(), p.data(), aux.data());
    plogp = cblas_ddot(GetNumberOfStrings(), p.data(), 1, aux.data(), 1);

    q.resize(GetNumberOfStrings());
    logq.resize(GetNumberOfStrings());

    if (HasUniquePaths())
        relative_path_probs.clear();
    else
        relative_path_probs.resize(GetNumberOfPaths());

    M.Init(GetNumberOfStrings(), GetNumberOfPaths(), Mrow.data(), Mcol.data(), ones.data());
    P.Init(GetNumberOfPaths(), GetNumberOfParameters(), Prow.data(), Pcol.data(), Pdata.data());
    C.Init(GetNumberOfParameters(), GetNumberOfConstraints(), Crow.data(), Ccol.data(), ones.data());

    FinalizeCallback();
}

double Learner::LogAuxiliaryVolume() const
{
    return LogSimplexVolume(auxiliary_parameters);
}

double Learner::LogModelVolume() const
{
    return model_volume;
}

double Learner::LogVolume() const
{
    return LogModelVolume() + LogAuxiliaryVolume();
}

double Learner::LogDetAuxiliaryHessian() const
{
    return aux_hessian;
}

MKL_INT Learner::GetNumberOfAuxParameters() const
{
    return (MKL_INT)auxiliary_parameters;
}

void Learner::ComputeModeledProbs()
{
    if (HasUniquePaths())
    {   // relative_path_probs holds nothing in this case!

        // logq = P.x
        P.dot(_x.data(), logq.data());

        // q = exp(P.x)
        vdExp(GetNumberOfStrings(), logq.data(), q.data());
    }else
    {
        aux.resize(GetNumberOfPaths());

        // relative_path_probs = P.x
        P.dot(_x.data(), relative_path_probs.data());

        // relative_path_probs = exp(P.x)
        vdExp(GetNumberOfPaths(), relative_path_probs.data(), relative_path_probs.data());

        // q = M.relative_path_probs
        M.dot(relative_path_probs.data(), q.data());

        // logq = log(q)
        vdLn(GetNumberOfStrings(), q.data(), logq.data());

        // aux <- M^t.q
        M.dot(q.data(), aux.data(), SPARSE_OPERATION_TRANSPOSE);

        // relative_path_probs /= M^t.q
        vdDiv(GetNumberOfPaths(), relative_path_probs.data(), aux.data(), relative_path_probs.data());
    }
}

void Learner::ComputeObjective()
{
    // p.log(p) - p.log(q)
    kl = plogp - cblas_ddot(GetNumberOfStrings(), p.data(), 1, logq.data(), 1);
}

MKL_INT Learner::GetNumberOfStrings() const
{
    return (MKL_INT)p.size();
}

MKL_INT Learner::GetNumberOfPaths() const { return (MKL_INT)Mcol.size(); }
MKL_INT Learner::GetNumberOfParameters() const { return (MKL_INT)Ccol.size(); }
MKL_INT Learner::GetNumberOfConstraints() const { return Ccol.empty() ? (MKL_INT)0 : Ccol.back() + 1; }

double Learner::GetCommonSupport() const
{
    return common_support;
}

double Learner::gKLDistance() const
{
    return GetKLDistance() + mxlogx(common_support);
}

void Learner::Init(int flags, const double* initialx)
{
    if (initialx)
        std::copy(initialx, initialx + GetNumberOfParameters(), _x.begin());
    InitCallback(flags);
}

void Learner::InitCallback(int)
{
}

bool Learner::HasUniquePaths() const
{
    return Mrow.size() == Mcol.size() + 1;
}

double Learner::GetKLDistance() const
{
    return kl;
}

const double* Learner::GetWeights() const
{
    return _x.data();
}
