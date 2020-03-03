#include "Transducer.h"

#include <string>
#include <limits>
#include <unordered_map>
#include <iostream>

#include "Utils.h"

Transducer::Transducer() { }

Transducer::Transducer(std::istream & is)
{
    Read(is);
}

static auto& mycast = SaturateCast<TransducerIndex>::template Do<size_t>;

static inline bool str_ends(TransducerIndex x)
{
    for (size_t i = 0; i < sizeof(TransducerIndex); ++i)
    {
        if (((const char*)(&x))[i] == '\0')
            return true;
    }
    return false;
}

static size_t str_space_required(const char * s)
{
    auto original = s;
    while (*s)
    {
        ++s;
    }
    return Round<sizeof(TransducerIndex)>::Do((s - original) + 1);
}

static void append_str(const char * s, std::vector<TransducerIndex>& v)
{
    const auto end = v.size();
    v.resize(end + str_space_required(s) / sizeof(TransducerIndex), 0);
    char* target = (char*)(&(v[end]));
    while (*s)
    {
        *target++ = *s++;
    }
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
        first_line.push_back(c);
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

    if (fread(&start_state_start, sizeof(start_state_start), 1, f) != 1)
        throw MyError("Invalid binary file format!");
    if (fread(&start_state_end, sizeof(start_state_end), 1, f) != 1)
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
    
    if (fwrite(&start_state_start, sizeof(start_state_start), 1, f) != 1)
        throw MyError("Cannot Write binary file!");
    if (fwrite(&start_state_end, sizeof(start_state_end), 1, f) != 1)
        throw MyError("Cannot Write binary file!");
    
    if (fwrite(transitions_table.data(), sizeof(TransducerIndex), transitions_table.size(), f) != transitions_table.size())
        throw MyError("Cannot Write binary file!");
}

void Transducer::Read(std::istream & is)
{
    std::unordered_map<TransducerIndex, std::array<TransducerIndex, 2>> state_pointers;

    transitions_table.clear();
    n_transitions = 0;
    n_states = 0;

    TransducerIndex previous_state = std::numeric_limits<TransducerIndex>::max();
    std::string line;

    Counter<std::string, TransducerIndex> states;

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
        std::string input;
        switch (parts.size())
        {
        case 1:
        case 2:
            // final state
            to = std::numeric_limits<TransducerIndex>::max();
            break;
        case 4:
        case 5:
            // ordinary transition
            to = states[parts[1]];
            input = parts[2];
            if (input == "@0@" || input == "@_EPSILON_SYMBOL_@")
                // these will end up deleting from input tape anyway
                // we don't care what will happen to output tape!
                input.clear();
            break;
        };
        if (previous_state != from)
        {
            state_pointers[from][0] = mycast(transitions_table.size());
            state_pointers[previous_state][1] = state_pointers[from][0];
            previous_state = from;
        }
        transitions_table.emplace_back(mycast(n_transitions));
        transitions_table.emplace_back(from);
        transitions_table.emplace_back(to);

        // input
        append_str(input.c_str(), transitions_table);

        ++n_transitions;
    }
    n_states = state_pointers.size();
    state_pointers[previous_state][1] = mycast(transitions_table.size());

    for (auto i = transitions_table.begin(); i < transitions_table.end(); ++i)
    {
        const auto& to = *(i + 2);
        if (to == std::numeric_limits<TransducerIndex>::max())
        {//final state has a MAXINT destination
            *(i + 1) = std::numeric_limits<TransducerIndex>::max();
        }else
        {
            *(i + 1) = state_pointers[to][0];
            *(i + 2) = state_pointers[to][1];
        }
        
        i += 3;
        while (!str_ends(*i))
        {
            ++i;
        }
    }

    start_state_start = state_pointers[0][0];
    start_state_end = state_pointers[0][1];
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

void Transducer::Lookup(const char* s, ResultHandler resulth, double time_l, size_t max_r, bool utf8)
{
    max_results = max_r;
    resulthandler = resulth;
    time_limit = time_l;
    clock.Tick();
    n_results = 0;
    path.clear();

    if (utf8)
        GetNextChar = GetNextUtf8Character;
    else
        GetNextChar = Next;

    lookup(s, start_state_start, start_state_end);
}

// ISSUE implementations don't follow HFST's own definition!
// @[PNDRCU][.][A-Z]+([.][A-Z]+)?@

void Transducer::lookup(const char* s, TransducerIndex beg, const TransducerIndex end)
{
    for (; beg < end; ++beg)
    {   // try outgoing edges
        const auto id = transitions_table[beg];
        const auto to_beg = transitions_table[beg + 1];
        const auto to_end = transitions_table[beg + 2];
        const auto input = reinterpret_cast<const char*>(transitions_table.data() + beg + 3);
        
        if (to_end == std::numeric_limits<TransducerIndex>::max() ||
            to_beg == std::numeric_limits<TransducerIndex>::max())
        {   //final state
            if (*s == '\0')
            {   // that's a result
                resulthandler(path);
            }
        }
        else if (StrEq2(input, "@_UNKNOWN_SYMBOL_@") || StrEq2(input, "@_IDENTITY_SYMBOL_@"))
        {   // consume one character
            path.emplace_back(id);
            lookup(GetNextUtf8Character(s), to_beg, to_end);
        }
        // 
        // epsilon is handled with a simple empty string
        // 
        else if (input[0] == '@' && input[2] == '.' && 
            (input[1] == 'P' || input[1] == 'N' || input[1] == 'D' || input[1] == 'R' || input[1] == 'C' || input[1] == 'U') &&
            input[3] != '\0')
        {   // TODO implement flags
            path.emplace_back(id);
            lookup(s, to_beg, to_end);
        }
        else if (const char* next = ContainsPrefix2(s, input))
        {   // a lead to follow
            path.emplace_back(id);
            lookup(next, to_beg, to_end);
        }

        beg += 3;
        while (!str_ends(transitions_table[beg]))
        {
            ++beg;
        }
    }

    if (!path.empty())
        path.pop_back();
}
