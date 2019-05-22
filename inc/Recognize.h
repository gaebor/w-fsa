#pragma once

#include <queue>
#include <functional>

#include "Utils.h"
#include "Fsa.h"

template<class Path>
class Recognizer
{
public:
    typedef std::function<Path&(Path&, const Fsa::Transitions::value_type&, const Fsa::Emissions::value_type&)> Accumulator;
    typedef std::function<void(const Path&)> ResultHandler;
private:
    Accumulator accumulator;
    ResultHandler resultHandler;

    const char* end_state;    

    struct QueueType
    {
        const char* word; //! what's left of the word
        const Fsa::State* state; //! current state in the automaton
        Path history; //! accumulated path
        QueueType(const char* w = "", const Fsa::State* s = nullptr, const Path& h = Path())
            : word(w), state(s), history(h)
        {
        }
    };
public:
    Recognizer(const char* endstate, Accumulator accumulator, ResultHandler resultHandler)
    :   accumulator(accumulator), resultHandler(resultHandler), end_state(endstate) {}

    void RecognizeDFS(const char* word, const Fsa::State& state, Path history = Path())const
    {
        for (const auto& transition : state.transitions)
        {
            if (StrEq::streq(transition.next->first, end_state))
            {   // reached end
                if (word[0] == '\0')
                {   //the string has been consumed
                    auto path(history);
                    accumulator(path, transition, Fsa::Emissions::value_type());
                    resultHandler(path);
                }
                continue;
            }
            const auto& nextstates = transition.next->second;
            for (const auto& nextemission : nextstates.emissions)
            {
                if (ContainsPrefix(word, nextemission.str))
                {   // a lead to follow
                    auto newpath(history);
                    accumulator(newpath, transition, nextemission);
                    RecognizeDFS(word + strlen(nextemission.str), nextstates, newpath);
                }
            }
        }
    }

    void RecognizeBFS(const char* original_word, const Fsa::State& state, Path history = Path())const
    {
        std::queue<QueueType> queue;
        queue.emplace(original_word, &state, history);

        while (!queue.empty())
        {
            auto& working = queue.front();
            auto current_word = working.word;
            for (const auto& transition : working.state->transitions)
            {
                if (StrEq::streq(transition.next->first, end_state))
                {   // reached end
                    if (current_word[0] == '\0')
                    {   //the string has been consumed
                        auto path(working.history);
                        accumulator(path, transition, Fsa::Emissions::value_type());
                        resultHandler(path);
                    }
                    continue;
                }
                const auto& nextstates = transition.next->second;
                for (const auto& nextemission : nextstates.emissions)
                {
                    if (ContainsPrefix(current_word, nextemission.str))
                    {   // a lead to follow
                        queue.emplace(current_word + strlen(nextemission.str), &nextstates, working.history);
                        accumulator(queue.back().history, transition, nextemission);
                    }
                }
            }
            queue.pop();
        }
    }

    template<bool BFS=true>
    void Recognize(const char* word, const Fsa::State& state, Path history = Path())const
    {
        if (BFS)
            RecognizeIndexBFS(word, state, history);
        else
            RecognizeIndexDFS(word, state, history);
    }
};
