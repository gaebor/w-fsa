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
    : common_support(0.0), plogp(0.0), kl(0.0),
    unique_path(false)
{
}

//Learner::Learner(const Fsa& fsa)
//:   common_support(0.0), plogp(0.0), kl(0.0),
//    unique_path(false)
//{
//    BuildConstraints(fsa); 
//}
//
//Learner::Learner(const std::string& filename)
//{
//    LoadMatrices(filename);
//}

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
        return false;
    ifs.open(filename + ".M");
    if (ifs)
    {
        ReadCsrMtx(ifs, ones, Mrow, Mcol);
        ifs.close();
    }
    else
        return false;
    ifs.open(filename + ".P");
    if (ifs)
    {
        ReadCsrMtx(ifs, Pdata, Prow, Pcol);
        ifs.close();
    }
    else
        return false;
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
            return false;
        ifs.close();
    }
    else
        return false;

    if (p.size() + 1 != Mrow.size())
        throw LearnerError("Size mismatch: ", p.size() + 1, " != ", Mrow.size());

    if (*std::max_element(Pcol.begin(), Pcol.end()) >= GetNumberOfParameters())
        throw LearnerError("Size mismatch: more path indexes than parameters in the automaton!");

    if (Mcol.back() + 1 != Prow.size() - 1)
        throw LearnerError("M cols != P rows");

    _x.assign(GetNumberOfParameters(), 0.0);
    unique_path = GetNumberOfStrings() == GetNumberOfPaths();

    return true;
}

std::vector<double> Learner::GetOptimizationInfo()
{
    return std::vector<double>();
}

std::string Learner::GetOptimizationHeader()
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

            Crow.emplace_back(Ccol.size());
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
    Crow.emplace_back(Ccol.size());
}

typedef std::vector<std::pair<MKL_INT, double>> Path;

void Learner::BuildFrom(const Fsa& fsa, const Corpus& corpus, bool bfs)
{
    BuildConstraints(fsa);
    BuildPaths(fsa, corpus, bfs);
}

void Learner::BuildPaths(const Fsa& fsa, const Corpus& corpus, bool bfs)
{
    auxiliary_parameters = 0;
    const auto& startState = fsa.GetTransitionMtx().at(fsa.GetStartState());
    bool has_path = false;
    Corpus::const_iterator wordIt; // currently processed word

    Recognizer<Path>::Accumulator accumulator(
    [endstate=fsa.GetEndState()](Path& history, const Fsa::Transitions::value_type& transition, const Fsa::Emissions::value_type& emission) -> Path&
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
            Mrow.emplace_back(Mcol.size());
            common_support += wordIt->second;
            p.emplace_back(wordIt->second);
        }
        has_path = true;

        Prow.emplace_back(Pcol.size());
        for (const auto& variable : path)
        {
            Pdata.emplace_back(variable.second);
            Pcol.emplace_back(variable.first);
        }

        Mcol.emplace_back(Mcol.size());        
    });
    Recognizer<Path> recognizer(fsa.GetEndState(), accumulator, resultHandler);
    auto recognize = bfs ? &Recognizer<Path>::RecognizeBFS : &Recognizer<Path>::RecognizeDFS;

    Pdata.clear(); Pcol.clear(); Prow.clear();
    Mcol.clear(); Mrow.clear();
 
    common_support = 0.0;
    p.clear();
    p.reserve(corpus.size()); // just an estimate

    size_t current_word = 0;
    ProgressIndicator(current_word, &current_word, 100.0 / corpus.size(),
        "\rRecognize: %4.3g%% ",
        [&]()
    {
        aux_hessian = 0;
        for (wordIt = corpus.begin(); wordIt != corpus.end(); ++wordIt)
        {
            ++current_word;
            has_path = false;
            (recognizer.*recognize)(wordIt->first.c_str(), startState, Path());
            
            if (!has_path)
            {   // word was not recognized
                ++auxiliary_parameters;
                aux_hessian -= std::log(wordIt->second);
            }
        }
    });

    if (Mcol.size() == Mrow.size())
        unique_path = true;

    Mrow.emplace_back(Mcol.size());
    Prow.emplace_back(Pcol.size());
}

void Learner::FinalizeCallback() {}

void Learner::Finalize()
{
    aux.resize(GetNumberOfStrings());

    ones.assign(std::max(GetNumberOfPaths(), GetNumberOfParameters()), 1.0);

    vdLn(GetNumberOfStrings(), p.data(), aux.data());
    plogp = cblas_ddot(p.size(), p.data(), 1, aux.data(), 1);

    q.resize(GetNumberOfStrings());
    logq.resize(GetNumberOfStrings());

    if (!unique_path)
        relative_path_probs.resize(GetNumberOfPaths());
    else
        relative_path_probs.clear();

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

size_t Learner::GetNumberOfAuxParameters() const
{
    return auxiliary_parameters;
}

void Learner::ComputeModeledProbs()
{
    if (unique_path)
    {   // relative_path_probs holds nothing in this case!

        // logq = P.x
        P.dot(_x.data(), logq.data());

        // q = exp(P.x)
        vdExp(logq.size(), logq.data(), q.data());
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
        vdLn(q.size(), q.data(), logq.data());

        // aux <- M^t.q
        M.dot(q.data(), aux.data(), SPARSE_OPERATION_TRANSPOSE);

        // relative_path_probs /= M^t.q
        vdDiv(GetNumberOfPaths(), relative_path_probs.data(), aux.data(), relative_path_probs.data());
    }
}

void Learner::ComputeObjective()
{
    // q.1
    // modeled_prob = cblas_ddot(GetNumberOfStrings(), q.data(), 1, ones.data(), 1);
    // p.log(p) - p.log(q)
    kl = plogp - cblas_ddot(GetNumberOfStrings(), p.data(), 1, logq.data(), 1);
}

size_t Learner::GetNumberOfStrings() const
{
    return p.size();
}

size_t Learner::GetNumberOfPaths() const { return Mcol.size(); }
size_t Learner::GetNumberOfParameters() const { return Ccol.size(); }
size_t Learner::GetNumberOfConstraints() const { return Ccol.empty() ? 0 : Ccol.back() + 1; }

double Learner::GetCommonSupport() const
{
    return common_support;
}

//double Learner::GetTotalModeledProb() const
//{
//    return modeled_prob;
//}

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

bool Learner::HasUniquePaths() const { return unique_path; }

double Learner::GetKLDistance() const
{
    return kl;
}

const double* Learner::GetWeights() const
{
    return _x.data();
}
