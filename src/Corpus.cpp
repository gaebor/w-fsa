#include "Corpus.h"

#include <unordered_set>

Corpus::Corpus()
{
}

void Corpus::Read(FILE * input)
{
    std::string word;
    clear();
    std::unordered_set<std::string> words;
    std::vector<char> content;

    if (!ReadContent(input, content))
        throw CorpusError("Cannot read file!");
    char* c = content.data();
    std::pair<const char*, char> result;

    result = GetWord(c, "\n");
    separator = result.first;
    if (separator.empty())
        separator = " ";

    ProgressIndicator(content.data(), &c, 100.0/content.size(), "\rCorpus: %4.3g%% ",
        [&]()
    {
        while (result.second != '\0')
        {
            word.clear();
            bool empty = true;
            do
            {
                result = GetWord(c, separator.c_str());
                if (result.second == '\n' || result.second == '\0')
                {
                    if (!empty)
                    {
                        if (words.insert(word).second)
                        {
                            emplace_back(word, atof(result.first));
                        }
                        else
                            throw CorpusError("\"", word, "\" is duplicate!");
                    }
                    break;
                }
                empty = false;
                word += result.first;
            } while (result.second);
        }
    });
    for (auto& wordIt : *this)
    {
        if (!std::isnormal(wordIt.second) || wordIt.second < 0)
        {
            throw CorpusError("\"", wordIt.first, "\" has probability ", wordIt.second, "!");
        }
    }
}

Corpus::~Corpus()
{
}

void Corpus::Renormalize()
{
    const double s = Sum();
    for (auto& word : *this)
        word.second /= s;
}

double Corpus::Sum() const
{
    double s = 0.0;
    for (const auto& word : *this)
        s += word.second;
    return s;
}
