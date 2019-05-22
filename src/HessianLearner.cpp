#include "HessianLearner.h"

#include "mkl.h"

#define _USE_MATH_DEFINES 
#include <math.h>

#include <unordered_set>
#include <algorithm>
#include <iterator>

HessianLearner::HessianLearner()
    : Learner(), solver_handle(0), include_Hf(false)
{
}

HessianLearner::HessianLearner(const Fsa & fsa)
:   Learner(fsa), solver_handle(nullptr), include_Hf(false)
{
}

HessianLearner::HessianLearner(const std::string & filename)
    : Learner(filename), solver_handle(nullptr), include_Hf(false)
{
}

HessianLearner::HessianLearner(const HessianLearner & other)
    : Learner(other), solver_handle(0)
{
    rhs = other.rhs;
    H = other.H;
    expx = other.expx;
    jgl = other.jgl;
    delta_x = other.delta_x;

    Hrow = other.Hrow; Hcol = other.Hcol;
    perm = other.perm;

    solver_opt = other.solver_opt;
    include_Hf = other.include_Hf;
}

void HessianLearner::BuildPathsCallback()
{
    rhs.resize(GetNumberOfAugmentedParameters());
    _x.resize(GetNumberOfAugmentedParameters(), 1.0);
    expx.resize(GetNumberOfParameters());
    jgl.resize(GetNumberOfParameters());
    delta_x.resize(GetNumberOfAugmentedParameters());
}

void HessianLearner::InitDss(bool reorder)
{
    solver_opt = MKL_DSS_ZERO_BASED_INDEXING +
        MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR +
        MKL_DSS_REFINEMENT_ON;
    
    MKL_INT error;
    if ((error = dss_create(solver_handle, solver_opt)) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Failed to create DSS! Error code: ", error);
    }
    const MKL_INT rows = GetNumberOfAugmentedParameters();
    const MKL_INT cols = GetNumberOfAugmentedParameters(), nnz = H.size();
    solver_opt = MKL_DSS_SYMMETRIC;
    if ((error = dss_define_structure(solver_handle, solver_opt, Hrow.data(), rows, cols, Hcol.data(), nnz)) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Unable to define structure! Error code: ", error);
    }
    solver_opt = reorder ? MKL_DSS_GET_ORDER : MKL_DSS_MY_ORDER;
    perm.resize(rows);
    if (!reorder)
    {   // stick to original order
        for (MKL_INT i = 0; i < rows; ++i)
            perm[i] = i;
    }
    if ((error = dss_reorder(solver_handle, solver_opt, perm.data())) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Unable to find reorder! Error code: ", error);
    }
}

HessianLearner::~HessianLearner()
{
    DeleteSolver();
}

void HessianLearner::OptimizationStep(double eta)
{
    ComputeRhs();
    ComputeObjective();
    
    // H = 0
    cblas_dscal(H.size(), 0.0, H.data(), 1);

    if (!unique_path && include_Hf)
        ComputeHf();

    //TODO determinant of H_f

    ComputeHg();
    
    // sol = H \ rhs
    MKL_INT nRhs = 1, error;
    solver_opt = MKL_DSS_INDEFINITE;
    if ((error = dss_factor_real(solver_handle, solver_opt, H.data())) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Unable to factor coefficient matrix! Error code: ", error);
    }

    solver_opt = MKL_DSS_REFINEMENT_ON;
    if ((error = dss_solve_real(solver_handle, solver_opt, rhs.data(), nRhs, delta_x.data())) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Unable to solve equation! Error code: ", error);
    }
    cblas_daxpy(delta_x.size(), -eta, delta_x.data(), 1, _x.data(), 1);
}

void HessianLearner::GetInertia(std::ostream & os)
{
    MKL_INT error;
    solver_opt = MKL_DSS_DEFAULTS;
    double ret_values[3];
    if ((error = dss_statistics(solver_handle, solver_opt, "Inertia", ret_values)) != MKL_DSS_SUCCESS)
    {
        throw LearnerError("Unable to get inertia! Error code: ", error);
    }
    os << '+' << size_t(ret_values[0]) << " -" << size_t(ret_values[1]) << " 0" << size_t(ret_values[2]) << std::endl;
}

void HessianLearner::InitCallback(int flags)
{
    size_t i = 0;
    size_t place = 4;
    ProgressIndicator<size_t>(0, &i, 1,
        "\rInitialize %5X ",
        [&]()
    {
        const size_t done = 0xE;
        const size_t skipped = 0xF;
        
        if (flags & 1)
        {
            _x.assign(GetNumberOfParameters(), 0.0);
            _x.resize(GetNumberOfParameters() + GetNumberOfConstraints(), 1.0);
            i += (done << (4 * place));
        }else
            i += (skipped << (4 * place));
        --place;
        
        if (flags & 2)
        {
            Renormalize();
            i += (done << (4 * place));
        }else
            i += (skipped << (4 * place));
        --place;
        
        if (flags & 4)
        {
            InitSlackVariables();
            i += (done << (4 * place));
        }else
            i += (skipped << (4 * place));
        --place;
        
        include_Hf = flags & 8;
        AssembleH(include_Hf);
        i += ((include_Hf ? done : skipped) << (4 * place));
        --place;
        
        InitDss(flags & 16);
        i += (((flags & 16) ? done : skipped) << (4 * place));
    });
}

size_t HessianLearner::GetNumberOfAugmentedParameters() const
{
    return GetNumberOfParameters() + GetNumberOfConstraints();
}

double HessianLearner::DeNormalizedFactor() const
{
    const auto n = GetNumberOfParameters(), k = GetNumberOfConstraints();
    return std::abs(rhs[n + cblas_idamax(k, rhs.data() + n, 1)]);
}

double HessianLearner::GradientError() const
{
    return std::abs(rhs[cblas_idamax(GetNumberOfParameters(), rhs.data(), 1)]);
}

size_t HessianLearner::GetNnzHessian() const
{
    return Hcol.size();
}

double HessianLearner::GetHessianFillRatio() const
{
    return (double)Hcol.size() / (Hrow.size() - 1);
}

double HessianLearner::ComputeLogDetHessian()
{
    // Hf is not zero, but hasn't been included
    if (!include_Hf && !unique_path)
        AssembleH(true);

    const auto n = GetNumberOfParameters();
    std::vector<MKL_INT> newHcol; newHcol.reserve(n); // at least
     // modify the nnz structure of H, you have to just erase the last element from each row
    for (size_t i = 0; i < n; ++i)
    {
        newHcol.insert(newHcol.end(), &Hcol[Hrow[i]], &Hcol[Hrow[i + 1] - 1]);
        Hrow[i] -= i;
    }
    Hrow.resize(n + 1);
    Hrow.back() = Hrow[n - 1] + 1;
    std::swap(newHcol, Hcol);
    H.resize(Hrow.back());

    // H = 0
    cblas_dscal(H.size(), 0, H.data(), 1);
    
    // fill in the values, similar to Hf

    ComputeExpX(); // these are the actual variables of the current Hessian

    static std::vector<std::pair<MKL_INT, double>> grad_qi;
    static std::vector<std::pair<std::pair<MKL_INT, MKL_INT>, double>> Hqijk;
    
    for (size_t str_idx = 0; str_idx < Mrow.size() - 1; ++str_idx)
    {
        const auto npaths = Mrow[str_idx + 1] - Mrow[str_idx];
        if (npaths == 1)
        {   // only diagonal elements
            const auto path_idx = Mrow[str_idx];
            for (auto rj = Prow[path_idx]; rj < Prow[path_idx + 1]; ++rj)
            {   // variables in the one and only path
                // p_i * #(x_j in path of p_i) / x_j^2
                
                GetCoord(Hrow, Hcol, H, Pcol[rj], Pcol[rj]) += p[str_idx] * Pdata[rj] / (expx[Pcol[rj]] * expx[Pcol[rj]]);
            }
        }
        else
        {   // equivocal str
            const auto piqi = p[str_idx] / q[str_idx];

            grad_qi.clear();
            for (auto path_idx = Mrow[str_idx]; path_idx < Mrow[str_idx + 1]; ++path_idx)
            {   // path
                const auto path_prob = path_probs[path_idx];

                Hqijk.clear();
                for (auto rj = Prow[path_idx]; rj < Prow[path_idx + 1]; ++rj)
                {   // variables in the path

                    // #(x_j in path of p_i) / x_j^2
                    SortedInsert(Hqijk, std::make_pair(Pcol[rj], Pcol[rj])) += Pdata[rj] / (expx[Pcol[rj]] * expx[Pcol[rj]]);

                    // path_prob * #{x_j in path} / exp(x_j)
                    SortedInsert(grad_qi, Pcol[rj]) += (path_prob / expx[Pcol[rj]]) * Pdata[rj];

                    for (auto rk = rj; rk < Prow[path_idx + 1]; ++rk)
                    {
                        // #{x_j in path}/exp(x_j) * #{x_k in path}/exp(x_k)
                        SortedInsert(Hqijk, std::make_pair(Pcol[rj], Pcol[rk])) -= (Pdata[rj] / expx[Pcol[rj]]) * (Pdata[rk] / expx[Pcol[rk]]);
                    }
                }
                cblas_dscal(Hqijk.size(), piqi * path_prob, &Hqijk.front().second, 3);
                for (const auto& jk : Hqijk)
                {
                    const auto j = jk.first.first;
                    const auto k_ = jk.first.second;
                    GetCoord(Hrow, Hcol, H, j, k_) += jk.second;
                }
            }
            cblas_dscal(grad_qi.size(), sqrt(p[str_idx]) / q[str_idx], &grad_qi.front().second, 2);
            for (size_t rj = 0; rj < grad_qi.size(); ++rj)
            {
                const auto j = grad_qi[rj].first;
                for (size_t rk = rj; rk < grad_qi.size(); ++rk)
                {
                    const auto k_ = grad_qi[rk].first;
                    GetCoord(Hrow, Hcol, H, j, k_) += grad_qi[rj].second * grad_qi[rk].second;
                }
            }
        }
    }
    return CalculateLogDetH(false);
}

void HessianLearner::PrintH(std::ostream & os) const
{
    PrintCsrMtx(os, H, Hrow, Hcol);
}

void HessianLearner::PrintEq(std::ostream & os) const
{
    PrintCsrMtx(os, H, Hrow, Hcol, rhs);
}

void HessianLearner::AssembleH(bool Hf)
{
    std::set<std::pair<MKL_INT, MKL_INT>> indices;
    std::vector<MKL_INT> equivocal_indices;

    // gather nnz indexes
    if (Hf && !unique_path)
    {
        for (size_t str_idx = 0; str_idx < Mrow.size() - 1; ++str_idx)
        {
            if (Mrow[str_idx] + 1 < Mrow[str_idx + 1])
            {   // equivocal str
                equivocal_indices.clear();
                auto first_path_idx = Mrow[str_idx];
                auto end_path_idx = Mrow[str_idx + 1];
                auto first_path_first_index = Prow[first_path_idx];
                auto end_path_end_index = Prow[end_path_idx];

                for (auto j = Pcol.begin() + first_path_first_index; j != Pcol.begin() + end_path_end_index; ++j)
                    SortedInsert(equivocal_indices, *j);
                
                for (auto i = equivocal_indices.begin(); i != equivocal_indices.end(); ++i)
                {
                    for (auto j = i; j != equivocal_indices.end(); ++j)
                        indices.emplace(*i, *j);
                }
            }
        }
    }

    // H_g
    for (size_t i = 0; i < GetNumberOfParameters(); ++i)
    {
        indices.emplace(i, i);
        indices.emplace(i, GetNumberOfParameters() + Ccol[i]);
    }
    for (size_t i = GetNumberOfParameters(); i < GetNumberOfAugmentedParameters(); ++i)
        indices.emplace(i, i);

    // convert to csr format
    H.clear(); Hrow.clear(); Hcol.clear();

    MKL_INT row = -1;
    for (const auto& ij : indices)
    {
        if (ij.first > row)
        {
            Hrow.emplace_back(Hcol.size());
            row = ij.first;
        }
        Hcol.emplace_back(ij.second);
        H.emplace_back(0.0);
    }
    Hrow.emplace_back(Hcol.size());
}

void HessianLearner::ComputeHf()
{
    static std::vector<std::pair<MKL_INT, double>> grad_qi;
    static std::vector<std::pair<std::pair<MKL_INT, MKL_INT>, double>> Hqijk;

    for (size_t str_idx = 0; str_idx < Mrow.size() - 1; ++str_idx)
    {
        const auto npaths = Mrow[str_idx + 1] - Mrow[str_idx];
        if (npaths > 1)
        {   // equivocal str
            const auto piqi = p[str_idx] / q[str_idx];

            grad_qi.clear();
            for (auto path_idx = Mrow[str_idx]; path_idx < Mrow[str_idx + 1]; ++path_idx)
            {   // path
                const auto path_prob = path_probs[path_idx];

                Hqijk.clear();
                for (auto rj = Prow[path_idx]; rj < Prow[path_idx + 1]; ++rj)
                {   // variables in the path
                    
                    // #{x_j in path} * path_prob
                    SortedInsert(grad_qi, Pcol[rj]) += path_prob * Pdata[rj];
                    
                    for (auto rk = rj; rk < Prow[path_idx + 1]; ++rk)
                    {
                        // #{x_j in path} * #{x_k in path}
                        SortedInsert(Hqijk, std::make_pair(Pcol[rj], Pcol[rk])) += Pdata[rj] * Pdata[rk];
                    }
                }
                static_assert(sizeof(decltype(Hqijk)::value_type) == 3 * sizeof(decltype(grad_qi)::value_type::second_type), "sizeof(pair<pair<MKL_INT, MKL_INT>, double>) != 3 * sizeof(double)!");
                cblas_dscal(Hqijk.size(), piqi * path_prob, &Hqijk.front().second, 3);
                for (const auto& jk : Hqijk)
                {
                    const auto j = jk.first.first;
                    const auto k_ = jk.first.second;
                    // H[j,k] -= pi / q_i * path_prob * Sum( #{x_j in path} * #{x_k in path} )
                    GetCoord(Hrow, Hcol, H, j, k_) -= jk.second;
                }
            }
            static_assert(sizeof(decltype(grad_qi)::value_type) == 2 * sizeof(decltype(grad_qi)::value_type::second_type), "sizeof(pair<MKL_INT, double>) != 2 * sizeof(double)!");
            // sqrt(qi)/qi
            cblas_dscal(grad_qi.size(), sqrt(p[str_idx]) / q[str_idx], &grad_qi.front().second, 2);
            for (size_t rj = 0; rj < grad_qi.size(); ++rj)
            {
                const auto j = grad_qi[rj].first;
                for (size_t rk = rj; rk < grad_qi.size(); ++rk)
                {
                    const auto k_ = grad_qi[rk].first;
                    // H[j,k] += pi * grad_qi_j/qi * grad_qi_k/qi
                    GetCoord(Hrow, Hcol, H, j, k_) += grad_qi[rj].second * grad_qi[rk].second;
                }
            }
        }
    }
}

void HessianLearner::ComputeExpX()
{
    vdExp(GetNumberOfParameters(), _x.data(), expx.data());
}

void HessianLearner::ComputeJgL()
{
    const MKL_INT n = GetNumberOfParameters(), k = GetNumberOfConstraints();
    double alpha = 1.0, beta = 0.0;
    mkl_dcsrmv("n", &n, &k, &alpha, "GxxCx",
        expx.data(), Ccol.data(), Crow.data(), Crow.data() + 1,
        _x.data() + n, &beta, jgl.data());
}

void HessianLearner::InitSlackVariables()
{
    double alpha = 1.0, beta = 0.0;
    const MKL_INT n = GetNumberOfParameters(), k = GetNumberOfConstraints();
    std::vector<double> sumexp2x(k);
    // this is going to contain the pseudo inverse of J
    auto& PInvJ = jgl;

    ComputeExpX();
    ComputeGrad();
    
    {   //PInvJ <- exp(2*x)

        // PInvJ <- _x
        cblas_dcopy(n, _x.data(), 1, PInvJ.data(), 1);
        // PInvJ *= 2
        cblas_dscal(n, 2.0, PInvJ.data(), 1);
        // PInvJ <- exp(PInvJ)
        vdExp(n, PInvJ.data(), PInvJ.data());
    }
    {
        // sumexp2x <- C^t.PInvJ
        mkl_dcsrmv("t", &n, &k, &alpha, "GxxCx",
            ones.data(), Ccol.data(), Crow.data(), Crow.data() + 1,
            PInvJ.data(), &beta, sumexp2x.data());
        // PInvJ <- C*sumexp2x
        mkl_dcsrmv("n", &n, &k, &alpha, "GxxCx",
            ones.data(), Ccol.data(), Crow.data(), Crow.data() + 1,
            sumexp2x.data(), &beta, PInvJ.data());
    }
    // PInvJ <- exp(x) / PInvJ
    vdDiv(n, expx.data(), PInvJ.data(), PInvJ.data());

    // -PInvJ.gradf
    alpha = -1.0; beta = 0.0;
    mkl_dcsrmv("t", &n, &k, &alpha, "GxxCx",
        PInvJ.data(), Ccol.data(), Crow.data(), Crow.data() + 1,
        rhs.data(), &beta, _x.data() + n);
}

void HessianLearner::ComputeGrad()
{
    static std::vector<double> aux, aux2;
    if (aux.empty() && unique_path) // grad_f is constant, we calculate it for the first time only
    {
        // aux <- -P^t p
        aux.resize(GetNumberOfParameters());
        const MKL_INT row = GetNumberOfPaths(); // equal the number of strings
        const MKL_INT col = GetNumberOfParameters();
        const double alpha = -1.0, beta = 0.0;
        mkl_dcsrmv("t", &row, &col, &alpha, "GxxCx",
            Pdata.data(), Pcol.data(), Prow.data(), Prow.data() + 1,
            p.data(), &beta, aux.data());
    }
    else
        aux.resize(GetNumberOfPaths()); // for future use

    ComputeModeledProbs();

    if (unique_path)
    {
        // rhs <- aux
        cblas_dcopy(GetNumberOfParameters(), aux.data(), 1, rhs.data(), 1);
    }
    else
    {
        MKL_INT row, col;
        aux2.resize(GetNumberOfStrings());

        // aux2 <- p/q
        vdDiv(GetNumberOfStrings(), p.data(), q.data(), aux2.data());

        // aux <- (diag(path_probs).M^t).p/q
        row = GetNumberOfStrings();
        double alpha = 1.0, beta = 0.0;
        col = GetNumberOfPaths();
        mkl_dcsrmv("t", &row, &col, &alpha, "GxxCx",
            path_probs.data(), Mcol.data(), Mrow.data(), Mrow.data() + 1,
            aux2.data(), &beta, aux.data());

        // rhs <- -P^t.aux
        alpha = -1.0;  beta = 0.0;
        row = GetNumberOfPaths();
        col = GetNumberOfParameters();
        mkl_dcsrmv("t", &row, &col, &alpha, "GxxCx",
            Pdata.data(), Pcol.data(), Prow.data(), Prow.data() + 1,
            aux.data(), &beta, rhs.data());
    }
}

void HessianLearner::ComputeG()
{
    const MKL_INT n = GetNumberOfParameters();
    const MKL_INT k = GetNumberOfConstraints();
    double alpha = 1.0, beta = 0.0;
    // g <- C^t.exp(x)
    mkl_dcsrmv("t", &n, &k, &alpha, "GxxCx",
        ones.data(), Ccol.data(), Crow.data(), Crow.data() + 1,
        expx.data(), &beta, rhs.data() + n);
    // g -= 1
    cblas_daxpy(k, -1.0, ones.data(), 1, rhs.data() + n, 1);
}

void HessianLearner::DeleteSolver()
{
    if (solver_handle)
    {
        solver_opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR;
        auto error = dss_delete(solver_handle, solver_opt);
        if (error != MKL_DSS_SUCCESS)
            throw LearnerError("Cannot delete solver of equation, error: ", error);
        solver_handle = 0;
    }
}

double HessianLearner::CalculateLogDetH(bool reorder)
{
    if (unique_path)
    {   // product of diagonal
        double result = 0;
        for (size_t i = 0; i < GetNumberOfParameters(); ++i)
        {
            if (H[i] > 0.0)
                result += log(H[i]);
            else
                return atof("inf");
        }
        return result;
    }
    else
    {   // factorize, get det
        MKL_INT error;
        const MKL_INT n = Hrow.size() - 1;
        const MKL_INT nnz = Hcol.size();
        if (solver_handle)
            DeleteSolver();

        solver_opt = MKL_DSS_ZERO_BASED_INDEXING +
            MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR +
            MKL_DSS_REFINEMENT_ON;
        if ((error = dss_create(solver_handle, solver_opt)) != MKL_DSS_SUCCESS)
        {
            throw LearnerError("Failed to create DSS for det H! Error code: ", error);
        }

        solver_opt = MKL_DSS_SYMMETRIC;
        if ((error = dss_define_structure(solver_handle, solver_opt, Hrow.data(), n, n, Hcol.data(), nnz)) != MKL_DSS_SUCCESS)
        {
            throw LearnerError("Unable to define structure of log det Hessian! Error code: ", error);
        }
        solver_opt = reorder ? MKL_DSS_GET_ORDER : MKL_DSS_MY_ORDER;
        perm.resize(n);
        if (!reorder)
        {   // stick to original order
            for (MKL_INT i = 0; i < n; ++i)
                perm[i] = i;
        }
        if ((error = dss_reorder(solver_handle, solver_opt, perm.data())) != MKL_DSS_SUCCESS)
        {
            throw LearnerError("Unable to find reorder! Error code: ", error);
        }
        solver_opt = MKL_DSS_INDEFINITE;
        if ((error = dss_factor_real(solver_handle, solver_opt, H.data())) != MKL_DSS_SUCCESS)
        {
            throw LearnerError("Unable to factor log det Hessian! Error code: ", error);
        }
        solver_opt = MKL_DSS_DEFAULTS;
        double det[2];
        if ((error = dss_statistics(solver_handle, solver_opt, "Determinant", det)) != MKL_DSS_SUCCESS)
        {
            throw LearnerError("Unable to get determinant of Hessian! Error code: ", error);
        }
        else
        {
            const double& det_pow = det[0], &det_base = det[1];
            if (det_base <= 0.0)
                return atof("inf");
            return log(det_base) + det_pow * M_LN10;
        }
    }
}

void HessianLearner::ComputeRhs()
{
    ComputeExpX();
    ComputeJgL();
    ComputeGrad();
    ComputeG();

    cblas_daxpy(GetNumberOfParameters(), 1.0, jgl.data(), 1, rhs.data(), 1);
}

void HessianLearner::ComputeHg()
{
    const MKL_INT n = GetNumberOfParameters(), k = GetNumberOfConstraints();
    if (n > 0)
        H[0] += jgl[0];
    for (MKL_INT i = 0; i < n - 1; ++i)
    {
        const auto r = Hrow[i + 1];
        H[r - 1] += expx[i];
        H[r] += jgl[i + 1];
    }
    if (n > 0)
    {
        H[H.size() - k - 1] += expx.back();
    }
}
