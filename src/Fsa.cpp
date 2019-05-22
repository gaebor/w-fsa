#include "Fsa.h"

#include <numeric>
#include <algorithm>
#include <iterator>

Fsa::Fsa()
    : m1(0), m2(0), n(0),
    separator(""), start_state(""), end_state("")
{}

Fsa& Fsa::operator=(const Fsa& other)
{
    content = other.content;
    auto offset = content.empty() ? 0 : content.data() - other.content.data();
    // TODO copy everything 
    /*transition_probs.clear();
    for (const auto& state : other.transition_probs)
    {
        auto& this_state = transition_probs[state.first + offset];
        for (const auto& transition : state.second.transitions)
        {
            this_state.transitions.emplace_back(transition.);
        }
    }*/
    transition_probs = other.transition_probs;
    m1 = other.m1;
    m2 = other.m2;
    n = other.n;
    if (other.separator < other.content.data() + other.content.size() && other.separator >= other.content.data())
        separator = other.separator + offset;
    else
        separator = other.separator;
    start_state = other.start_state + offset;
    end_state = other.end_state + offset;
    
    return *this;
}

Fsa::Fsa(const Fsa & other)
{
    *this = other;
}

Fsa::~Fsa() {}

void Fsa::Dump(FILE * out) const
{
    fprintf(out, "%s\n%s\n%s\n", separator, start_state, end_state);

    for (const auto& state : transition_probs)
    {
        if (StrEq::streq(state.first, end_state))
            continue;
        fprintf(out, "%s", state.first);
        for (const auto& emission : state.second.emissions)
        {
            fprintf(out, "%s%s%s%g",
                separator, emission.str,
                separator, emission.logprob);
        }
        fprintf(out, "\n%s", state.first);
        for (const auto& transition : state.second.transitions)
        {
            fprintf(out, "%s%s%s%g",
                separator, transition.next->first,
                separator, transition.logprob);
        }
        fprintf(out, "\n");
    }
}

void Fsa::Read(FILE* input)
{
    Clear();

    if (!ReadContent(input, content))
        throw FsaError("Unable to read file!");
    const auto expected_size = AllocateStates();

    char* c = content.data();

    for (auto x : { &separator, &start_state, &end_state })
    {
        *x = GetWord(c, "\n").first;
    }
    if (IsEmpty(separator))
        separator = " ";

    for (auto x : { &start_state, &end_state })
    {
        if (ContainsPrefix(*x, separator))
            throw FsaError("Invalid FSA format! "
                "Start or end state contains the separator! ",
                "\"", separator, "\" is in \"", *x, "\"");
    }
    if (StrEq::streq(start_state, end_state))
        throw FsaError("Invalid FSA format! "
            "Start and end states should be different! ",
            "\"", start_state, "\"==\"", end_state, "\"");

    ProgressIndicator(c, &c, 100.0 / content.size(), "\rAutomaton: %4.3g%% ",
        [&]()
    {
        while (*c)
            ReadOneState(c);
    });
    if (transition_probs.size() > expected_size)
        throw FsaError("Invalid FSA format! "
            "There are more states than rows in the automaton file! ",
            transition_probs.size(), " > ", expected_size);
    AssignIndices();
}

size_t Fsa::AllocateStates()
{
    auto lines = std::accumulate(content.begin(), content.end(),
        size_t(0), [](size_t i, char c) { return i + (c == '\n' ? 1 : 0); });
    const auto expected_states = (std::max<size_t>(lines, 3) - 3) / 2 + 1;
    transition_probs.reserve(expected_states);
    transition_probs.max_load_factor(0.8f);
    return expected_states;
}

void Fsa::ReadOneState(char*& c)
{
    auto result = GetWord(c, separator);
    auto this_state = result.first;

    if (IsEmpty(this_state) || ContainsPrefix(this_state, end_state))
    {   // skip this line
        GetWord(c, "\n");
        return;
    }
    Keyed<double> emissions, transitions;
    // emissions
    do
    {
        result = GetWord(c, separator);
        auto word = result.first;
        if (StrEq::streq(this_state, start_state) && !IsEmpty(word))
            throw FsaError("Invalid FSA format! "
                "Start state should emit empty string instead of \"", word, "\"!");
        auto insertion = emissions.emplace(word, 0.0);
        if (!insertion.second)
            throw FsaError("Invalid FSA format! "
                "Emission \"", word, "\" of state \"", this_state, "\" appears more than once!");
        result = GetWord(c, separator);
        (insertion.first->second) = atof(result.first);
    } while (result.second != '\n' && result.second != '\0');

    // check for errors
    if (emissions.empty())
    {
        throw FsaError("Invalid FSA format! "
            "State \"", this_state, "\" should have positive number of emissions, even if empty emission!");
    }

    // transitions
    if (!(StrEq::streq(this_state, GetWord(c, separator).first)))
    {
        throw FsaError("Invalid FSA format! "
            "You should enlist transitions of \"", this_state,
            "\" after emissions of the same state!");
    }
    do
    {
        result = GetWord(c, separator);
        auto word = result.first;
        auto insertion = transitions.emplace(word, 0.0);
        if (!insertion.second)
            throw FsaError("Invalid FSA format! "
                "Transition \"", this_state, "\" -> \"", word,
                "\" appears more than once!");
        if (StrEq::streq(word, start_state))
            throw FsaError("Invalid FSA format! "
                "\"", this_state, "\" connects to start state \"",
                start_state, "\"!");
        result = GetWord(c, separator);
        (insertion.first->second) = atof(result.first);
    } while (result.second != '\n' && result.second != '\0');

    // check for errors
    if (transitions.empty())
    {
        throw FsaError("Invalid FSA format! "
            "State \"", this_state, "\" should have positive number of transitions!");
    }
    State new_state;
    for (const auto& emission : emissions)
    {
        new_state.emissions.emplace_back(
            emission.first,
            emission.second
            );
    }
    for (const auto& transition : transitions)
    {
        new_state.transitions.emplace_back(
            &(*transition_probs.emplace(transition.first, State()).first),
            transition.second
            );
    }
    transition_probs[this_state] = new_state;
}

void Fsa::AssignIndices()
{
    m1 = 0; m2 = 0; n = 0;
    
    for (auto& state : transition_probs)
    {
        auto& emissions = state.second.emissions;
        auto& transitions = state.second.transitions;
        
        if (emissions.size() == 1)
        {   // unequivocal
            emissions.begin()->index = -1;
        } else
        {   // equivocal
            for (auto& e : emissions)
                e.index = n++;
        }

        if (transitions.size() == 1)
        {   // unequivocal
            transitions.begin()->index = -1;
        }
        else
        {   // equivocal
            for (auto& t : transitions)
                t.index = n++;
        }

        m1 += transitions.size();
        m2 += emissions.size();
    }
}

size_t Fsa::GetNumberOfStates() const
{
    return transition_probs.size();
}

size_t Fsa::GetNumberOfTransitions() const
{
    return m1;
}

size_t Fsa::GetNumberOfEmissions() const
{
    return m2;
}

size_t Fsa::GetNumberOfFreeParameters() const
{
    return m1 + m2 - 2 * (transition_probs.size() - 1);
}

size_t Fsa::GetNumberOfParameters() const
{
    return n;
}

size_t Fsa::GetNumberOfConstraints() const
{
    return n - GetNumberOfFreeParameters();
}

void Fsa::Clear()
{
    m1 = 0; m2 = 0; n = 0;
    separator = ""; start_state = ""; end_state = "";
    transition_probs.clear();
}

const Fsa::TransitionMtx & Fsa::GetTransitionMtx() const
{
    return transition_probs;
}

void Fsa::ResetWeights(const double* x)
{
    for (auto& t : transition_probs)
    {
        auto& transitions = t.second.transitions;
        auto& emissions = t.second.emissions;

        for (auto& edge : emissions)
            edge.logprob = (edge.index >= 0 ? x[edge.index] : 0.0);

        for (auto& edge : transitions)
            edge.logprob = (edge.index >= 0 ? x[edge.index] : 0.0);
    }
}
//
//StrPaths Fsa::Recognize(const char* word)const
//{
//    StrPaths paths;
//    RecognizeDFS<StrPath>(word, transition_probs.at(start_state), StrPath(), paths,
//        [this](StrPath history, const Transitions::value_type& transition, const Emissions::value_type& emission)
//    {
//        history.emplace_back(emission.str, transition.next->first);
//        return history;
//    }
//    );
//    return paths;
//}
