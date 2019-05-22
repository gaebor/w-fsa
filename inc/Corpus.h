#pragma once

#include <cstdio>

#include <string>
#include <unordered_map>

class Corpus : public std::unordered_map<std::string, double>
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
