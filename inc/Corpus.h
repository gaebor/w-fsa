#pragma once

#include <cstdio>

#include <string>
#include <utility>
#include <vector>

#include "Utils.h"

struct CorpusError : public MyError
{
    using MyError::MyError;
};

class Corpus : public std::vector<std::pair<std::string, double>>
{
public:
    Corpus();
    void Read(FILE* input);
    ~Corpus();
    void Renormalize();
    double Sum()const;
private:
    std::string separator;
};
