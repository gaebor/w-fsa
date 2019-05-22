#pragma once

#include <vector>
#include <ostream>

#include "mkl_types.h"

#include "Fsa.h"
#include "Corpus.h"

struct LearnerError : public MyError
{
    using MyError::MyError;
};

class Learner
{
public:
    Learner(); 
    Learner(const Fsa& fsa);
    Learner(const std::string& filename);

    void BuildPaths(const Fsa& fsa, const Corpus& corpus, bool bfs = true);

    //!
    /*!
    builds C matrix, reads _x from FSA, calculates log-volume of the model

    +------+
    |1     |
    |1     |
    |1     |
    |  1   |
    |  1   |
    |    1 |
    |    1 |
    |    1 |
    +------+
    Rows: n (parameters)
    Cols: k (constraints)
    */
    void BuildConstraints(const Fsa& fsa);

    virtual ~Learner();

    //! makes sure weights sum up to one.
    void Renormalize();

    const double* GetWeights()const;

    void PrintC(std::ostream& os)const;
    void PrintJ(std::ostream& os)const;
    void PrintP(std::ostream& os)const;
    void PrintM(std::ostream& os)const;

    bool SaveMatrices(const std::string filename)const;

    bool LoadMatrices(const std::string filename);

    double GetCommonSupport()const;
    double GetTotalModeledProb()const;

    //! computes the modeled probabilities of paths, strings and log-prob of strings and stores them internally
    void ComputeModeledProbs();
    //! computes KL and sum(q)
    void ComputeObjective();

    size_t GetNumberOfStrings()const;
    size_t GetNumberOfPaths()const;
    size_t GetNumberOfParameters()const;
    size_t GetNumberOfConstraints()const;

    bool HasUniquePaths()const;

    double GetKLDistance()const;
    double gKLDistance()const;

    virtual void OptimizationStep(double eta = 1.0) = 0;
    virtual void Init(int flags, const double* initialx = nullptr);

    double LogModelVolume()const;
    double LogAuxiliaryVolume()const;
    double LogVolume()const;
    double LogDetAuxiliaryHessian()const;
    size_t GetNumberOfAuxParameters()const;
protected:
    virtual void BuildPathsCallback();
    virtual void InitCallback(int flags);
private:
    /*!
    sets vector p and matrices P and M
    computes common_support
    calculates number of auxiliary parameters (unrecognized strings)
    calculates log det auxiliary Hessian
    finds out unique path property
    */
    void BuildStructureMtx(const Fsa& fsa, const Corpus& corpus, bool bfs = false);

    //!
    /*!
        computes plogp
        allocates memory to ones, q, log(q), path_probs
        assigns: ones
    */
    void AssignConstants();
protected:
    std::vector<double> _x; //!< parameters
    std::vector<double> p, q, logq, path_probs;
private:
    double common_support, plogp, kl, modeled_prob, aux_hessian, model_volume;
    size_t auxiliary_parameters;
protected:
    std::vector<double> ones;
    //! constraints matrix
    /*!
        number of model parameters rows, number of constraints columns
        C[i, j] = 1.0 if x_i appears in the j-th constraint
    */
    std::vector<MKL_INT> Crow, Ccol;
    //! paths matrix
    /*!
        number of paths rows, number of model parameters columns
        P[i, j] = how many times the x_j parameter participated in the i-th path
    */
    std::vector<MKL_INT> Prow, Pcol;
    std::vector<double> Pdata;

    //! 
    /*!
        number of paths cols, number of string rows
        This matrix is identity <-> unique path property
    */
    std::vector<MKL_INT> Mrow, Mcol;
protected:
    bool unique_path;
};
