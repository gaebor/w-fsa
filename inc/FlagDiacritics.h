#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <type_traits>

class FlagDiacritics
{
public:
    enum OpType { Pop, Nop, Rop, Dop, Cop, Uop };
    // which feature, which value, + or -, 0 if value is empty
    typedef std::unordered_map<std::string, std::pair<std::string, int>> State;
    typedef std::tuple<OpType, std::string, std::string> Operation;
    static bool IsIt(const char* input);
    static Operation Parse(const char* flags);
    void Read(const char* flags);
    // bool Apply(const char* flags);

    // if false, then state is not modified,
    // if true, then operation is applied to the state 
    static bool Apply(const Operation& op, State& state);
    static bool Apply(const char* flags, State& state);

protected:
    std::pair<unsigned char, unsigned char> GetFeatureValue(const char* feature, const char* value)const
    {
        std::pair<unsigned char, unsigned char> p;
        for (auto& f : flag_map)
        {
            if (f.first == feature)
            {
                p.second = 1;
                for (auto& v : f.second)
                {
                    if (v == value)
                        return p;
                    ++p.second;
                }
            }
            ++p.first;
        }
        return p;
    }
private:
    std::unordered_map<std::string, std::unordered_set<std::string>> flag_map;
    State flag_state;
};
