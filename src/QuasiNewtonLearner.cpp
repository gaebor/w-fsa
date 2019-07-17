#include "QuasiNewtonLearner.h"

#include "mkl.h"
#include "Utils.h"

#include <math.h>

QuasiNewtonLearner::QuasiNewtonLearner()
    : Learner(), fixed_lambda(false)
{
}

QuasiNewtonLearner::~QuasiNewtonLearner()
{
}

void QuasiNewtonLearner::FinalizeCallback()
{
    grad.resize(GetNumberOfParameters());
    expx.resize(GetNumberOfParameters());
    rhs.resize(GetNumberOfParameters());

    lambda.assign(GetNumberOfConstraints(), 1.0);
    g.resize(GetNumberOfConstraints());

    Jg.Init(GetNumberOfParameters(), GetNumberOfConstraints(), Crow.data(), Ccol.data(), expx.data());
}

void QuasiNewtonLearner::InitCallback(int flags)
{
    if (flags & 1)
    {
        _x.assign(GetNumberOfParameters(), 0.0);
    }
    
    if (flags & 2)
    {
        Renormalize();
    }
    
    fixed_lambda = (flags & 4) != 0;
}

void QuasiNewtonLearner::ComputeExpX()
{
    vdExp(GetNumberOfParameters(), _x.data(), expx.data());
}

std::string QuasiNewtonLearner::GetOptimizationHeader()const
{
    return   "       KL"
            "   graderr"
            "     g_min"
            "     g_max"
            " lambdamin"
            "      rmin";
}

std::vector<double> QuasiNewtonLearner::GetOptimizationInfo()
{
    std::vector<double> result(7);
    result[0] = GetKLDistance();

    result[1] = grad_error;

    result[2] = g_min;
    result[3] = g_max;
    
    result[4] = lambda_min;
    
    if (!HasUniquePaths())
    {
        result[6] = (double)cblas_idamin(GetNumberOfPaths(), relative_path_probs.data(), 1);
        result[5] = relative_path_probs[(MKL_INT)result[6]];
    }
    return result;
}

bool QuasiNewtonLearner::HaltCondition(double tol)
{
    return grad_error <= tol && std::abs(g_min) <= tol && std::abs(g_max) <= tol;
}

void QuasiNewtonLearner::ComputeGrad()
{
    if (grad_aux.empty())
    {   // initialize
        if (HasUniquePaths())
        {   // gradient is constant: -P^t.p
            grad_aux.resize(GetNumberOfParameters());
            P.dot(p.data(), grad_aux.data(), SPARSE_OPERATION_TRANSPOSE, -1.0, 0.0);
        }
        else
        {   // a useful quantity is precomputed: -M^t.p
            grad_aux.resize(GetNumberOfPaths());
            M.dot(p.data(), grad_aux.data(), SPARSE_OPERATION_TRANSPOSE, -1.0, 0.0);
        }
    }

    ComputeModeledProbs();
    if (HasUniquePaths())
    {
        // grad <- grad_aux
        cblas_dcopy(GetNumberOfParameters(), grad_aux.data(), 1, grad.data(), 1);
    }
    else
    {
        aux.resize(GetNumberOfPaths());

        // aux <- (-M^t.p)*relative_path_probs
        vdMul(GetNumberOfPaths(), relative_path_probs.data(), grad_aux.data(), aux.data());

        // grad <- P^t.aux
        P.dot(aux.data(), grad.data(), SPARSE_OPERATION_TRANSPOSE);
    }
}

void QuasiNewtonLearner::ComputeLambdaNext()
{
    lambda_min = *std::min_element(lambda.begin(), lambda.end());

    if (fixed_lambda)
    {
        C.dot(grad.data(), lambda.data(), SPARSE_OPERATION_TRANSPOSE, -1.0, 0.0);
    }
    else
    {// lambda <- 1/(g+1)*(g*lambda -C.grad)
        // lambda <- g*lambda
        vdMul(GetNumberOfConstraints(), lambda.data(), g.data(), lambda.data());
        // lambda -= C.grad
        C.dot(grad.data(), lambda.data(), SPARSE_OPERATION_TRANSPOSE, -1.0, 1.0);
        // g <- g+1
        vdAdd(GetNumberOfConstraints(), g.data(), ones.data(), g.data());

        vdDiv(GetNumberOfConstraints(), lambda.data(), g.data(), lambda.data());
    }
}

void QuasiNewtonLearner::ComputeG()
{
    const MKL_INT k = GetNumberOfConstraints();

    // g = -1
    g.assign(k, -1.0);

    // g += C^t.exp(x)
    C.dot(expx.data(), g.data(), SPARSE_OPERATION_TRANSPOSE, 1.0, 1.0);

    g_min = *std::min_element(g.begin(), g.end());
    g_max = *std::max_element(g.begin(), g.end());
}

void QuasiNewtonLearner::OptimizationStep(double eta, bool )
{
    ComputeExpX();
    ComputeG();
    ComputeGrad();
    ComputeObjective();

    aux.resize(GetNumberOfParameters());

    {// rhs = grad + J_g.lambda
        // aux <- J_g.lambda
        Jg.dot(lambda.data(), aux.data());

        // rhs <- grad
        cblas_dcopy(GetNumberOfParameters(), grad.data(), 1, rhs.data(), 1);
        // rhs += aux
        cblas_daxpy(GetNumberOfParameters(), 1.0, aux.data(), 1, rhs.data(), 1);
    }
    grad_error = std::abs(rhs[cblas_idamax(GetNumberOfParameters(), rhs.data(), 1)]);

    ComputeLambdaNext();

    {
        // rhs <- grad
        cblas_dcopy(GetNumberOfParameters(), grad.data(), 1, rhs.data(), 1);

        // rhs += J_g.lambda_next
        Jg.dot(lambda.data(), rhs.data(), SPARSE_OPERATION_NON_TRANSPOSE, 1.0, 1.0);
    }

    // rhs <- rhs/A
    vdDiv(GetNumberOfParameters(), rhs.data(), aux.data(), rhs.data());
    
    // x += -eta*rhs
    cblas_daxpy(GetNumberOfParameters(), -eta, rhs.data(), 1, _x.data(), 1);

}
