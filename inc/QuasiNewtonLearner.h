#pragma once

#include "Learner.h"

#include "mkl_types.h"

/*!
the following will be explicitly solved, without intel MKL
<pre>
              H              _x   rhs
+---------------+--------+  +--+  +--+
|               |        |  |  |  |  |
|               |        |  |  |  |  |
|   H_g         |   J_g  |  |x |  |grad_f + J_g*l
|               |        |  |  |  |  |
|               |        |* |  | =|  |
+---------------+--------+  |--|  |--|
|               |        |  |  |  |  |
|    J_g^T      |   0    |  |l |  |g(x)
|               |        |  |  |  |  |
+---------------+--------+  +--+  +--+
</pre>
*/
class QuasiNewtonLearner : public Learner
{
public:
    QuasiNewtonLearner();

    virtual ~QuasiNewtonLearner();

    virtual void OptimizationStep(double eta = 1.0, bool verbose = false);

    virtual std::vector<double> GetOptimizationInfo();
    virtual std::string GetOptimizationHeader()const;
    virtual bool HaltCondition(double tol);
protected:
    virtual void FinalizeCallback();
    virtual void InitCallback(int flags);
    
    void ComputeExpX();

    //! same as HessianLearner::ComputeGrad
    /*! 
        calls ComputeModeledProbs
    */
    void ComputeG();
    void ComputeGrad();
    //! lambda <- lambda_next
    /*! requires g, grad
        fucks up g in the meantime
    */
    void ComputeLambdaNext();
private:

    std::vector<double> grad, expx, lambda, g, grad_aux, rhs;
    SparseMtxHandle Jg;

    double grad_error, lambda_min, g_min, g_max;
    bool fixed_lambda;
};
