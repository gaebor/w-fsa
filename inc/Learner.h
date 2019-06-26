#pragma once

#include <vector>
#include <ostream>

#include "mkl_types.h"

#include "Fsa.h"
#include "Corpus.h"
#include "Utils.h"

struct LearnerError : public MyError
{
    using MyError::MyError;
};

class Learner
{
public:
    Learner(); 
    //Learner(const Fsa& fsa);
    //Learner(const std::string& filename);

    void BuildFrom(const Fsa& fsa, const Corpus& corpus, bool bfs = true);

    virtual ~Learner();

    //! makes sure weights sum up to one.
    void Renormalize();

    void RewriteWeights(Fsa& fsa)const;
    const double* GetWeights()const;


    void PrintC(FILE* f)const;
    void PrintJ(FILE* f)const;
    void PrintP(FILE* f)const;
    void PrintM(FILE* f)const;

    bool SaveMatrices(const std::string& filename)const;

    bool LoadMatrices(const std::string& filename);

    virtual std::vector<double> GetOptimizationInfo();
    virtual std::string GetOptimizationHeader();
    virtual std::vector<double> GetOptimizationResult(bool verbose= false);

    virtual bool HaltCondition(double tol);

    double GetCommonSupport()const;
    // double GetTotalModeledProb()const;

    //! computes the modeled (log-)probabilities and the relative_path_probs if needed (!unique_path)
    void ComputeModeledProbs();
    //! computes KL and sum(q)
    void ComputeObjective();

    MKL_INT GetNumberOfStrings()const;
    MKL_INT GetNumberOfPaths()const;
    MKL_INT GetNumberOfParameters()const;
    MKL_INT GetNumberOfConstraints()const;

    bool HasUniquePaths()const;
     
    double GetKLDistance()const;
    double gKLDistance()const;

    virtual void OptimizationStep(double eta = 1.0, bool verbose = false) = 0;
    virtual void Init(int flags, const double* initialx = nullptr);

    double LogModelVolume()const;
    double LogAuxiliaryVolume()const;
    double LogVolume()const;
    double LogDetAuxiliaryHessian()const;
    MKL_INT GetNumberOfAuxParameters()const;
    
    //!
    /*!
        computes plogp
        allocates memory to ones, q, log(q), path_probs
        assigns: ones
    */
    void Finalize();
protected:
    virtual void FinalizeCallback();
    virtual void InitCallback(int flags);
protected:
    //! parameters
    std::vector<double> _x; 
    //! measured probabilities 
    std::vector<double> p;
    //! modeled probabilities 
    std::vector<double> q, logq;
    /*! if unique_path, then relative_path_probs holds some auxiliary data!
        otherwise holds the relative path probabilities
    */
    std::vector<double> relative_path_probs;
    //! for various purposes, but only used locally
    std::vector<double> aux;
private:
    double common_support, plogp, kl, aux_hessian, model_volume;
    size_t auxiliary_parameters;

    /*!
    sets vector p and matrices P and M
    computes common_support
    calculates number of auxiliary parameters (unrecognized strings)
    calculates log det auxiliary Hessian
    finds out unique path property
    */
    void BuildPaths(const Fsa& fsa, const Corpus& corpus, bool bfs = true);

    /*! Trims unnecessary variables from P and C matrices
        it should be run after BuildPaths
    */
    void Trim();

    double GetWeight(MKL_INT i)const;

    std::vector<MKL_INT> trimmed_weights;

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
protected:
    std::vector<double> ones;
    //! constraints matrix
    /*!
        number of model parameters rows, number of constraints columns
        C[i, j] = 1.0 if x_i appears in the j-th constraint
    */
    std::vector<MKL_INT> Crow, Ccol;
    SparseMtxHandle C;

    //! paths matrix
    /*!
        number of paths rows, number of model parameters columns
        P[i, j] = how many times the x_j parameter participated in the i-th path
    */
    std::vector<MKL_INT> Prow, Pcol;
    std::vector<double> Pdata;
    SparseMtxHandle P;
    //! 
    /*!
        number of paths cols, number of string rows
        This matrix is identity <-> unique path property
    */
    std::vector<MKL_INT> Mrow, Mcol;
    SparseMtxHandle M;

    std::vector<std::pair<MKL_INT, std::vector<MKL_INT>>> equivocal_str_indices;
};
