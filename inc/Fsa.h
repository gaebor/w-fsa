#pragma once

#include <cstdio>

#include <vector>
#include <list>
#include <functional>
#include <set>

#include "Utils.h"

#include "mkl_types.h"

#ifndef FSA_CHECK_STATES
#define FSA_CHECK_STATES false
#endif

struct FsaError : public MyError
{
    using MyError::MyError;
};

class Fsa
{
public:
    struct NamedProb
    {
        const char* str;
        double logprob;
        MKL_INT index;
        NamedProb(const char* s = "", double l = 0.0, MKL_INT i = -1)
        : str(s), logprob(l), index(i)
        {}
    };
    struct IndexedProb
    {
        double logprob;
        MKL_INT index;
        IndexedProb(double l = 0.0, MKL_INT i = -1)
        : logprob(l), index(i)
        {}
    };
    struct State;
    typedef Keyed<State> TransitionMtx;

    // https://github.com/gaebor/CppPurgatory/blob/master/src/class_inside_class.cpp
    // https://stackoverflow.com/questions/41173415/self-referential-use-of-unordered-map-causes-problems-for-gcc-5-3-but-not-clang
    // error: 'std::pair<_T1, _T2>::second' has incomplete type
    // typedef const TransitionMtx::value_type* NextPtr;
    typedef const std::pair<const CStr, State>* NextPtr;
    struct NextState
    {
        NextPtr next;
        double logprob;
        MKL_INT index;
        NextState(NextPtr n, double l = 0.0, MKL_INT i = -1)
        : next(n), logprob(l), index(i)
        {}
    };
    typedef std::vector<NamedProb> Emissions;
    typedef std::vector<NextState> Transitions;
    struct State
    {
        Emissions emissions;
        Transitions transitions;
    };
public:
    Fsa();
    Fsa(const Fsa& other);
    Fsa& operator=(const Fsa& other);
    ~Fsa();

    void Dump(FILE* out)const;

    void Read(FILE* input);

    size_t GetNumberOfStates()const;
    size_t GetNumberOfTransitions()const;
    size_t GetNumberOfEmissions()const;
    size_t GetNumberOfFreeParameters()const;
    size_t GetNumberOfParameters()const;
    size_t GetNumberOfConstraints()const;

    const char* GetStartState()const { return start_state; }
    const char* GetEndState()const { return end_state; }

    const TransitionMtx& GetTransitionMtx()const;
    TransitionMtx& GetTransitionMtx();

private:
    void AssignIndices();
    //! reads how many states are expected and reserves so that no rehash will happen
    size_t AllocateStates();
    void Clear();
    void ReadOneState(char*& content);

    TransitionMtx transition_probs;

    size_t m1; //!< number of transition edges
    size_t m2; //!< number of emission edges
    size_t n; //!< number of free parameters

    CStr separator, start_state, end_state;

    std::vector<char> content;
};
