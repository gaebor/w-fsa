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

#ifdef _MSC_VER
#   define _ftell _ftelli64
#   define _fseek _fseeki64
#else
#   define _ftell ftello
#   define _fseek fseeko
#endif

void Transducer::ReadBinary(FILE* f)
{
    transitions_table.clear();
    n_transitions = 0;
    n_states = 0;

    std::string first_line;
    int c;
    while (( c = getc(f)) != EOF && c != 0)
    {
        first_line.push_back(static_cast<char>(c));
    }
    std::istringstream iss(first_line);
    size_t s;
    iss >> s;
    if (s != sizeof(TransducerIndex))
    {
        throw MyError("Binary file does not have the same alignment as compiled code (", s, "!=", sizeof(TransducerIndex), ")!");
    }
    if (c == EOF)
        return;

    if (fread(start_state.data(), sizeof(start_state), 1, f) != 1)
        throw MyError("Invalid binary file format!");
    s = _ftell(f);
    if (_fseek(f, 0, SEEK_END))
        throw MyError("Seeking failed during reading binary file!");
    transitions_table.resize((_ftell(f) - s) / sizeof(TransducerIndex));
    if (_fseek(f, s, SEEK_SET))
        throw MyError("Seeking failed during reading binary file!");
    if (fread(transitions_table.data(), sizeof(TransducerIndex), transitions_table.size(), f) != transitions_table.size())
        throw MyError("Cannot read transition table from binary file!");
}

void Transducer::DumpBinary(FILE* f)
{
    std::ostringstream oss;
    oss << sizeof(TransducerIndex);
    if (fwrite(oss.str().c_str(), 1, oss.str().size() + 1, f) != oss.str().size() + 1)
        throw MyError("Cannot Write binary file!");
    
    if (fwrite(start_state.data(), sizeof(start_state), 1, f) != 1)
        throw MyError("Cannot Write binary file!");
    
    if (fwrite(transitions_table.data(), sizeof(TransducerIndex), transitions_table.size(), f) != transitions_table.size())
        throw MyError("Cannot Write binary file!");
}

void Transducer::Read(std::istream & is)
{
    std::unordered_map<TransducerIndex, ToPointers> state_pointers;

    transitions_table.clear();
    n_transitions = 0;
    n_states = 0;

    TransducerIndex previous_state = std::numeric_limits<TransducerIndex>::max();
    std::string line;

    Counter<std::string, TransducerIndex> states;
    {

    while (std::getline(is, line))
    {
        std::vector<std::string> parts;
        size_t pos = 0, pos_next;
        while (std::string::npos != (pos_next = line.find('\t', pos)))
        {
            parts.emplace_back(line.begin() + pos, line.begin() + pos_next);
            pos = pos_next + 1;
        }
        parts.emplace_back(line.begin() + pos, line.end());
        TransducerIndex from = states[parts[0]], to;
        std::string input, output;
        float weight = 0;
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
    n_results = 0;
    path.clear();
    myclock.Tick();
    flag_failed = false;

    switch (strategy)
    {
    case FlagStrategy::IGNORE:
        return lookup<IGNORE>(s, start_state[0], start_state[1]);
    case FlagStrategy::NEGATIVE:
        return lookup<NEGATIVE>(s, start_state[0], start_state[1]);
    default:
        return lookup<OBEY>(s, start_state[0], start_state[1]);
    };
}

}
