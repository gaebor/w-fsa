#include "atto/Transducer.h"

#include <string>
#include <limits>
#include <unordered_map>
#include <sstream>
#include <array>

#include "Utils.h"
#include "atto/FlagDiacritics.h"
#include "atto/Record.h"

namespace atto {

Transducer::Transducer() { }

Transducer::Transducer(std::istream & is)
{
    Read(is);
}

void Transducer::Read(std::istream & is)
{
    std::unordered_map<TransducerIndex, ToPointers> state_pointers;

    transitions_table.clear();
    n_transitions = 0;
    n_states = 0;

    TransducerIndex previous_state = 0;
    std::string line;
    std::vector<std::string> parts;
    Counter<std::string, TransducerIndex> states;
    std::string input, output;
    float weight;

    while (std::getline(is, line))
    {
        parts.clear();
        size_t pos = 0, pos_next;
        while (std::string::npos != (pos_next = line.find('\t', pos)))
        {
            parts.emplace_back(line.begin() + pos, line.begin() + pos_next);
            pos = pos_next + 1;
        }
        parts.emplace_back(line.begin() + pos, line.end());
        TransducerIndex from = states[parts[0]], to;
        input.clear(); output.clear();
        weight = 0;
        switch (parts.size())
        {
        case 1:
        case 2:
            // final state
            to = std::numeric_limits<TransducerIndex>::max();
            if (parts.size() == 2)
                std::istringstream(parts[1]) >> weight;
            break;
        case 4:
        case 5:
            // ordinary transition
            to = states[parts[1]];
            input = parts[2];
            output = parts[3];
            for (auto s : {&input, &output})
            {
                if (*s == "@0@" || *s == "@_EPSILON_SYMBOL_@")
                    s->clear();
            }
            if (parts.size() == 5)
                std::istringstream(parts[4]) >> weight;
            break;
        default:
            throw MyError("AT&T file at line", n_transitions + 1, " has length ", parts.size());
        };
        if (previous_state != from)
        {
            state_pointers[from][0] = (TransducerIndex)transitions_table.size();
            state_pointers[previous_state][1] = state_pointers[from][0];
            previous_state = from;
        }
        
        //! write the binary format (pre-compile)
        transitions_table.emplace_back((TransducerIndex)n_transitions);
        transitions_table.emplace_back(from); // this will be to_start
        transitions_table.emplace_back(to); // this will be to_end
        transitions_table.emplace_back();
        static_assert(sizeof(weight) == sizeof(TransducerIndex), "");
        std::memcpy(&transitions_table.back(), &weight, sizeof(weight));

        if (FlagDiacritics<>::IsIt(input.c_str()))
        {
            if (input != output)
                    throw MyError("Input and output tape symbols doesn't math for flag diacritic: \"", input, "\" != \"", output, "\"");
            fd_table.Read(input.c_str());
        }
        CopyStr(input.c_str(), transitions_table);
        CopyStr(output.c_str(), transitions_table);
        
        ++n_transitions;
    }

    n_states = state_pointers.size();
    state_pointers[previous_state][1] = SaturateCast<TransducerIndex>::Do(transitions_table.size());

    fd_table.CalculateOffsets();

    // write the binary format
    const auto end = transitions_table.data() + transitions_table.size();
    for (RecordIterator<> i(transitions_table.data()); i < end; ++i)
    {
        auto& from = i.GetFrom();
        auto& to = i.GetTo();
     
        if (from != std::numeric_limits<TransducerIndex>::max() && 
            to != std::numeric_limits<TransducerIndex>::max())
        {   // non-final state
            from = state_pointers[to][0];
            to = state_pointers[to][1];
        }// these are 'pointers' now, not indexes
        
        if (fd_table.IsIt(i.GetInput()))
            fd_table.Compile(i.GetInput());
    }

    // supposing that the START is "0"
    start_state = state_pointers[states["0"]];
}

size_t Transducer::GetNumberOfTransitions() const
{
    return n_transitions;
}

size_t Transducer::GetNumberOfStates() const
{
    return n_states;
}

size_t Transducer::GetAllocatedMemory() const
{
    return sizeof(TransducerIndex)*transitions_table.size();
}

void Transducer::Lookup(const char* s, FlagStrategy strategy)
{
    auto lookup_fn = &Transducer::lookup<OBEY>;
    
    switch (strategy)
    {
    case FlagStrategy::IGNORE:
        lookup_fn = &Transducer::lookup<IGNORE>;
        break;
    case FlagStrategy::NEGATIVE:
        lookup_fn = &Transducer::lookup<NEGATIVE>;
        break;
    default:
        break;
    };

    (this->*lookup_fn)(nullptr, 0, 0);
    return (this->*lookup_fn)(s, start_state[0], start_state[1]);
}

}
