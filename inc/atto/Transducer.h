#pragma once

#include <vector>
#include <istream>
#include <functional>
#include <stdio.h>

#include "atto/FlagDiacritics.h"
#include "atto/char.h"
#include "Utils.h"

namespace atto {

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
    std::vector<TransducerIndex> transitions_table;
    TransducerIndex start_state_start, start_state_end;
    size_t n_transitions, n_states;
    
    template<FlagStrategy strategy>
    void lookup(const char* s, TransducerIndex beg, const TransducerIndex end)
    {
        if ((max_results == 0 || n_results < max_results) &&
            (max_depth == 0 || path.size() < max_depth) &&
            (time_limit == 0 || myclock.Tock() < time_limit))
            for (; beg < end; ++beg)
            {   // try outgoing edges
                const auto id = transitions_table[beg];
                const auto to_beg = transitions_table[beg + 1];
                const auto to_end = transitions_table[beg + 2];
                const auto input = reinterpret_cast<const char*>(transitions_table.data() + beg + 3);

                if (to_end == std::numeric_limits<TransducerIndex>::max() ||
                    to_beg == std::numeric_limits<TransducerIndex>::max())
                {   //final state
                    if (*s == '\0' && (strategy != NEGATIVE || flag_failed))
                    {   // that's a result
                        ++n_results;
                        path.emplace_back(id);
                        (*resulthandler)(path);
                        path.pop_back();
                    }
                }
                else if (StrEqual<Encoding::UTF8>(input, "@_UNKNOWN_SYMBOL_@") || 
                         StrEqual<Encoding::UTF8>(input, "@_IDENTITY_SYMBOL_@"))
                {   // consume one character
                    // TODO this can be hastened if the special symbols are shorter!
                    path.emplace_back(id);
                    lookup<strategy>(GetNextCharacter<Encoding::UTF8>(s), to_beg, to_end);
                }
                // 
                // epsilon is handled with a simple empty string
                // 
                else if (fd_table.IsIt(input))
                {
                    if (strategy == IGNORE)
                    {   // go with it, no matter what
                        path.emplace_back(id);
                        lookup<strategy>(s, to_beg, to_end);
                    }else
                    {
                        const auto flagsate = fd_state;
                        if (fd_table.Apply(input, fd_state))
                        {
                            path.emplace_back(id);
                            lookup<strategy>(s, to_beg, to_end);
                        }
                        else if (strategy == NEGATIVE)
                        {
                            const bool previous_fail = flag_failed;
                            flag_failed = true;
                            path.emplace_back(id);
                            lookup<strategy>(s, to_beg, to_end);
                            flag_failed = previous_fail;
                        }                            
                        fd_state = flagsate;
                    }
                }
                else if (const char* next = StrPrefix<Encoding::UTF8>(s, input))
                {   // a lead to follow
                    path.emplace_back(id);
                    lookup<strategy>(next, to_beg, to_end);
                }

                beg += 3;
                while (!StrEnds<UTF8>(transitions_table[beg]))
                {
                    ++beg;
                }
            }

        if (!path.empty())
            path.pop_back();
    }

    size_t n_results;
    Transducer::Path path;
    FlagDiacritics<>::State fd_state;
    FlagDiacritics<> fd_table;
    ::Clock<> myclock;
    bool flag_failed;
};

}
