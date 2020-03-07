#pragma once

#include <vector>
#include <istream>
#include <exception>
#include <functional>
#include <regex>
#include <unordered_map>
#include <array>
#include <stdio.h>

#include "Utils.h"
#include "FlagDiacritics.h"

#ifndef ATTOL_NUMBER
typedef unsigned int TransducerIndex;
#else
typedef ATTOL_NUMBER TransducerIndex;
#undef ATTOL_NUMBER
#endif

class Transducer
{
public:
    Transducer();
    Transducer(std::istream& is);
    void Read(std::istream& is);
    void ReadBinary(FILE* f);
    void DumpBinary(FILE* f);
    //! including to finishing from a final state
    size_t GetNumberOfTransitions()const;
    //! including start state
    size_t GetNumberOfStates()const;
    size_t GetAllocatedMemory()const;

    typedef std::vector<TransducerIndex> Path;
    typedef std::function<void(const Path&)> ResultHandler;

    void Lookup(const char* s, const ResultHandler& resulthandler, double time_limit = 0.0, size_t max_results = 0, size_t max_depth = 0);
private:
    std::vector<TransducerIndex> transitions_table;
    TransducerIndex start_state_start, start_state_end;
    size_t n_transitions, n_states;
    
    const ResultHandler* resulthandler;
    void lookup(const char* s, TransducerIndex beg, const TransducerIndex end);

};