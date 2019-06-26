#pragma once

#include "Learner.h"

#include "mkl_types.h"
#include "mkl_dss.h"

/*!
<pre>
              H              _x   rhs
+---------------+--------+  +--+  +--+
|               |        |  |  |  |  |
|               |        |  |  |  |  |
|   H_g + H_f   |   J_g  |  |x |  |grad_f + J_g*l
|               |        |  |  |  |  |
|               |        |* |  | =|  |
+---------------+--------+  |--|  |--|
|               |        |  |  |  |  |
|    J_g^T      |   0    |  |l |  |g(x)
|               |        |  |  |  |  |
+---------------+--------+  +--+  +--+
</pre>
*/
class HessianLearner : public Learner
{
public:
    HessianLearner();

    virtual ~HessianLearner();

    //! fills _x[n:]
    void InitSlackVariables();

    virtual void OptimizationStep(double eta = 1.0, bool verbose = false);

    MKL_INT GetNumberOfAugmentedParameters()const;

    //! returns how far the automaton is from a probabilistic one
    double DeNormalizedFactor()const;
    //! returns the max absolute value in the gradient vector
    double GradientError()const;

    MKL_INT GetNnzHessian()const;
    double GetHessianFillRatio()const;

    double ComputeLogDetHessian();

    void PrintH(FILE* f)const;
    void PrintEq(FILE* f)const;
    virtual std::vector<double> GetOptimizationInfo();
    virtual std::vector<double> GetOptimizationResult(bool verbose=false);
    virtual std::string GetOptimizationHeader();
    virtual bool HaltCondition(double tol);
protected:
    virtual void FinalizeCallback();
    virtual void InitCallback(int flags);

private:
    void InitDss(bool reorder=true);

    //!< structure of H
    void AssembleH(bool Hf = false);

    //! exp(x)
    void ComputeExpX();
    ////! J<sub>g</sub>.lambda
    //void ComputeJgL();

    //! rhs <- [grad f + J_g.lambda, g]
    /*!
        calls
            ComputeExpX
            ComputeGrad
            ComputeJgL
            ComputeG
    */
    void ComputeRhs();
    //! adds H_g to H
    /*! <i>H<sub>g</sub></i> is the following
    <pre>
        [------------------+-----]
        [                  |     ]
        [                  |     ]
        [ diag(J_g.lambda) |  J  ]
        [                  |     ]
        [                  |     ]
        [------------------+-----]
        [       J^t        |  0  ]
        [------------------+-----]
    </pre>
    */
    void ComputeHg();

    //! adds H<sub>f</sub> to the upper-left corner of H
    /*!
        does nothing if unique_path is true
    <pre>
                       ( grad q_i   (grad q_i)^t   H_q_i )
        H_f = Sum  p_i ( -------- . ------------ - ----- )
               i       (    q_i          q_i        q_i  )
       
        sum goes over equivocal strings

        P^t . diag(exp(P.x)) . M^t . diag(1/q) . diag(p) . diag(1/q) . M . diag(exp(P.x)) . P -
        - P^t . diag(P(str)*Q(path|str)) . P
                     |                |
                     +-M^t(p/q)*path--+
                       M^t(p/q)*(exp(P.x))
    </pre>
    */
    void ComputeHf();

    //! rhs[:n] <- grad f
    /*!
        calls ComputeModeledProbs
        
        M_ij = 1 if the j-th path participates in the i-th string, 0 otherwise
        </pre>
                       j   
            +---------------------+
            |1 1 1     |            |
            |     1 1  |            |
          i | - - -  1 1 1 1        |
            |                1      |
            |                  1 1 1|
            +-----------------------+

        P_jk = # {x_k occurred in the j-th path}
                      k
            +-------------+
            |         |   |
            |         |   |
            |         |   |
          j | - - - - #   |
            |             |
            |             |
            +-------------+

        path_prob = exp(P.x)
        q = M.exp(P.x)
        f = - p.log(q)
        grad f = Jac_q^t.(-p/q) =

            +--------------(grad q_i)/q_i-------+
            |                                   |
            +------- Jac^t ---------+           |
            |                       |           |
        = - P^t . diag(exp(P.x)). M^t . diag(1/q) . p
                  |   |        |                |
                  |   path probs                |
                  +------ relative path prob ---+
        </pre>
    */
    void ComputeGrad();

    //! rhs[n:] <- C^t.exp(x) - 1
    void ComputeG();

    std::vector<double> rhs, H, expx, grad_aux;
    SparseMtxHandle Jg;

    std::vector<MKL_INT> Hrow, Hcol;
    DssSolverHandler solver;
    MKL_INT solver_opt;
    double error, lambda_min;
    bool include_Hf, do_reorder, degenerate;
};
