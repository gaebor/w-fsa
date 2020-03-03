#include "FlagDiacritics.h"

#include "Utils.h"

bool FlagDiacritics::IsIt(const char * input)
{
    return (input[0] == '@' && input[2] == '.' &&
        (input[1] == 'P' || input[1] == 'N' || input[1] == 'D' || input[1] == 'R' || input[1] == 'C' || input[1] == 'U') &&
        input[3] != '\0');
}

bool FlagDiacritics::Apply(const Operation& op, State& flag_state)
{
    const auto& feature = std::get<1>(op);
    const auto& value = std::get<2>(op);
    switch (std::get<0>(op)) {
    case Pop: // positive set
        flag_state[feature].first = value;
        flag_state[feature].second = 1;
        return true;

    case Nop: // negative set
        flag_state[feature].second = -1;
        return true;

    case Rop: // require
        if (value.empty()) // empty require
            return !flag_state[feature].first.empty();
        else // nonempty require
            return (flag_state[feature].first == value && flag_state[feature].second == 1);

    case Dop: // disallow
        if (value.empty()) // empty disallow
            return flag_state[feature].first.empty();
        else // nonempty disallow
            return (flag_state[feature].first != value) || (flag_state[feature].second == -1);

    case Cop: // clear
        flag_state[feature].first.clear();
        flag_state[feature].second = 0;
        return true;

    case Uop: // unification
        if (flag_state[feature].first.empty()  || /* if the feature is unset or */
            (flag_state[feature].first == value && flag_state[feature].second == 1) || /* the feature is at
                                                  this value already
                                                  or */
            (flag_state[feature].second == -1 &&
            (flag_state[feature].first != value)) /* the feature is
                                                         negatively set
                                                         to something
                                                         else */
            )
        {
            flag_state[feature].first = value;
            flag_state[feature].second = 1;
            return true;
        }
        return false;
    }
    throw; // for the compiler's peace of mind
}

bool FlagDiacritics::Apply(const char * flags, State & state)
{
    return Apply(Parse(flags), state);
}

static FlagDiacritics::OpType char_to_operator(char c)
{
    switch (c) {
    case 'P': return FlagDiacritics::Pop;
    case 'N': return FlagDiacritics::Nop;
    case 'R': return FlagDiacritics::Rop;
    case 'D': return FlagDiacritics::Dop;
    case 'C': return FlagDiacritics::Cop;
    case 'U': return FlagDiacritics::Uop;
    default:
        throw MyError("Invalid flag diacritic operation: '", c, "'");
    }
}
void FlagDiacritics::Read(const char * flags)
{
    auto t = Parse(flags);
    if (std::get<2>(t).empty())
        flag_map[std::get<1>(t)];
    else
        flag_map[std::get<1>(t)].emplace(std::get<2>(t));
}

std::tuple<FlagDiacritics::OpType, std::string, std::string> 
    FlagDiacritics::Parse(const char * flags)
{
    std::string feature = flags + 3;
    std::string value;

    // @
    feature.pop_back();

    if (feature.find('.') != std::string::npos)
    {
        value = feature.substr(feature.find('.') + 1);
        feature = feature.substr(0, feature.find('.'));
    }
    return std::tuple<FlagDiacritics::OpType, std::string, std::string>(char_to_operator(flags[1]), feature, value);
    //if (value.empty())
    //    flag_map[feature];
    //else
    //    flag_map[feature].emplace(value);
}
