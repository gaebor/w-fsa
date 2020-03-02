#pragma once

#include <vector>
#include <istream>
#include <exception>
#include <functional>
#include <regex>
#include <unordered_map>
#include <array>

#include "Utils.h"

#ifndef ATTOL_NUMBER
typedef unsigned int TransducerIndex;
#else
typedef ATTOL_NUMBER TransducerIndex;
#undef ATTOL_NUMBER
#endif

size_t str_space_required(const char* s);

void append_str(const char* s, std::vector<TransducerIndex>& v);

bool str_ends(TransducerIndex x);

class Transducer
{
public:
    Transducer();
    Transducer(std::istream& is);
    void Read(std::istream& is);
    //! including to finishing from a final state
    size_t GetNumberOfTransitions()const;
    //! including start state
    size_t GetNumberOfStates()const;
    size_t GetAllocatedMemory()const;

    typedef std::vector<TransducerIndex> Path;
    typedef std::function<void(const Path&)> ResultHandler;

    void Lookup(const char* s, ResultHandler resulthandler, double time_limit = 0.0, size_t max_results = -1, bool utf8 = true);
private:
    std::vector<TransducerIndex> transitions_table;
    TransducerIndex start_state_start, start_state_end;
    TransducerIndex n_transitions, n_states;

    Path path;
    size_t max_results;
    size_t n_results;
    ResultHandler resulthandler;
    double time_limit;
    Clock<> clock;

    const char* (*GetNextChar)(const char*);

    static const char* Next(const char* x) { return x + 1; }

    void lookup(const char* s, TransducerIndex beg, const TransducerIndex end);

};