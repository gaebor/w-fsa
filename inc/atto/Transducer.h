#pragma once

#include <vector>
#include <istream>
#include <functional>
#include <array>
#include <tuple>
#include <stdio.h>

#include "atto/FlagDiacritics.h"
#include "atto/char.h"
#include "Utils.h"
#include "atto/Record.h"

namespace atto {

#ifndef ATTOL_NUMBER
typedef unsigned int TransducerIndex;
#else
typedef ATTOL_NUMBER TransducerIndex;
#undef ATTOL_NUMBER
#endif

class Transducer
{
private:
    struct ToPointers : std::array<TransducerIndex, 2>
    {
        ToPointers() { fill(0); }
    };

    std::vector<TransducerIndex> transitions_table;
    ToPointers start_state;
    size_t n_transitions, n_states;
    size_t n_results;
    
    FlagDiacritics<> fd_table;
    typedef typename decltype(fd_table)::State FlagState;
public:
    typedef std::vector<std::tuple<TransducerIndex, float, FlagState, const char*, const char*>> Path;
private:
    Transducer::Path path;
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

    typedef std::function<void(const Path&)> ResultHandler;

    size_t max_results;
    size_t max_depth;
    double time_limit;
    Transducer::ResultHandler* resulthandler;

    enum FlagStrategy
    {
        OBEY,
        IGNORE,
        NEGATIVE
    };

    void Lookup(const char* s, FlagStrategy strategy = FlagStrategy::OBEY);
private:
    template<FlagStrategy strategy>
    void lookup(const char* s, TransducerIndex beg, const TransducerIndex end)
    {
        if ((max_results == 0 || n_results < max_results) &&
            (max_depth == 0 || path.size() < max_depth) &&
            (time_limit == 0 || myclock.Tock() < time_limit))
        {
            for (RecordIterator<> i(transitions_table.data() + beg); i < transitions_table.data() + end; ++i)
            {   // try outgoing edges
                const auto input = i.GetInput();
                auto current_flag_state = path.empty() ? FlagState() : std::get<2>(path.back());

                if (i.GetTo() == std::numeric_limits<TransducerIndex>::max() ||
                    i.GetFrom() == std::numeric_limits<TransducerIndex>::max())
                {   //final state
                    if (*s == '\0' && (strategy != NEGATIVE || flag_failed))
                    {   // that's a result
                        ++n_results;
                        path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, input, i.GetOutput());
                        (*resulthandler)(path);
                        path.pop_back();
                    }
                }
                else if (StrEqual<Encoding::UTF8>(input, "@_UNKNOWN_SYMBOL_@") ||
                    StrEqual<Encoding::UTF8>(input, "@_IDENTITY_SYMBOL_@"))
                {   // consume one character
                    // TODO this can be hastened if the special symbols are shorter!
                    path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, input, i.GetOutput());
                    lookup<strategy>(GetNextCharacter<Encoding::UTF8>(s), i.GetFrom(), i.GetTo());
                }
                // 
                // epsilon is handled with a simple empty string
                // 
                else if (fd_table.IsIt(input))
                {
                    if (strategy == IGNORE)
                    {   // go with it, no matter what
                        path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, i.GetOutput(), i.GetOutput());
                        lookup<strategy>(s, i.GetFrom(), i.GetTo());
                    }
                    else
                    {
                        if (fd_table.Apply(input, current_flag_state))
                        {
                            path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, i.GetOutput(), i.GetOutput());
                            lookup<strategy>(s, i.GetFrom(), i.GetTo());
                        }
                        else if (strategy == NEGATIVE)
                        {
                            const bool previous_fail = flag_failed;
                            flag_failed = true;
                            path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, i.GetOutput(), i.GetOutput());
                            lookup<strategy>(s, i.GetFrom(), i.GetTo());
                            flag_failed = previous_fail;
                        }
                    }
                }
                else if (const char* next = StrPrefix<Encoding::UTF8>(s, input))
                {   // a lead to follow
                    path.emplace_back(i.GetId(), i.GetWeight(), current_flag_state, input, i.GetOutput());
                    lookup<strategy>(next, i.GetFrom(), i.GetTo());
                }
            }
        }
        if (!path.empty())
            path.pop_back();
    }
private:
    ::Clock<> myclock;
    bool flag_failed;
};

}
