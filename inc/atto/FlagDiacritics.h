#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <limits>
#include <type_traits>
#include <climits>

#include <Utils.h>

namespace atto {

template<class CharType = char, class FlagStorageType = int>
class FlagDiacritics
{
    static_assert(std::is_integral<CharType>::value && std::is_integral<FlagStorageType>::value,
        "atto::FlagDiacritics requires an integral type as template parameter!");
public:    
    typedef std::basic_string<CharType> string;
    typedef CharType* cstr;
    typedef const CharType* ccstr;

    struct State
    {
        State() : bitfield(0) {}
        //! returns a mask, where there are 1's in places [i, j)
        static inline FlagStorageType Mask(unsigned char i, unsigned char j)
        {
            // mask 0000001111000000
            //      5432109876543210
            //           j   i      
            return (((FlagStorageType)1 << (j - i)) - (FlagStorageType)1) << i;
        }
        FlagStorageType  Get(unsigned char i, unsigned char j)const
        {
            auto x = (((FlagStorageType)1 << (j - i)) - (FlagStorageType)1) & (bitfield >> i);
            if (x >> (j - i - 1))
            {   // negative value
                x |= (FlagStorageType)(-1) << (j - i);
            }
            return x;
        }
        void Set(unsigned char i, unsigned char j, FlagStorageType new_value)
        {
            const auto mask = Mask(i, j);
            // delete what's in place
            bitfield &= ~mask;
            // add new value to the right place
            bitfield |= mask & (new_value << i);
        }
        void clear()
        {
            bitfield = 0;
        }
    protected:
        FlagStorageType bitfield;
    };

    std::vector<FlagStorageType> GetValues(const State& s)const
    {
        std::vector<FlagStorageType> result;
        for (size_t i = 0; i + 1 < offsets.size(); ++i)
        {
            result.push_back(s.Get(offsets[i], offsets[i + 1]));
        }
        return result;
    }

    static bool IsIt(ccstr input)
    {
        return (input[0] == '@' && input[2] == '.' &&
            (input[1] == 'P' || input[1] == 'N' || input[1] == 'D' || input[1] == 'R' || input[1] == 'C' || input[1] == 'U') &&
            input[3] != 0);
    }
    static std::pair<string, string> Parse(ccstr flags)
    {
        std::pair<string, string> p;
        auto& feature = p.first;
        auto& value= p.second;
        feature = flags + 3;

        // @
        feature.pop_back();

        const size_t i = feature.find(CharType('.'));
        if (i != string::npos)
        {
            value = feature.substr(i + 1);
            feature = feature.substr(0, i);
        }
        return p;
    }

    void Read(ccstr flags)
    {
        const auto p = Parse(flags);
        if (p.second.empty())
            flag_map[p.first];
        else
            flag_map[p.first].emplace(p.second);
    }
    // replaces text encoded flag diacritic with binary optimized format
    void Compile(cstr flags)const
    {
        const auto p = Parse(flags);
        const auto b = GetFeatureValue(p.first.c_str(), p.second.c_str());
        flags[3] = b.first;
        flags[4] = b.second;
    }

    // if false, then state is not modified,
    // if true, then operation is applied to the state 
    bool Apply(ccstr flags, State& state)const
    {
        const unsigned char feature = flags[3];
        const unsigned char value = flags[4];
        const unsigned char bit_start = offsets[feature - 1];
        const unsigned char bit_end = offsets[feature];
        const auto current_value = state.Get(bit_start, bit_end);
        switch (flags[1])
        {
        case 'P': // positive set
            state.Set(bit_start, bit_end, value);
            return true;
        
        case 'N': // negative set
            state.Set(bit_start, bit_end, -value);
            return true;
        
        case 'R': // require
            if (value == 0) // empty require
                return current_value != 0;
            else // nonempty require
                return current_value == value;
        
        case 'D': // disallow
            if (value == 0) // empty disallow
                return current_value == 0;
            else // nonempty disallow
                return current_value != value;
        
        case 'C': // clear
            state.Set(bit_start, bit_end, 0);
            return true;
        
        case 'U': // unification
            if (current_value == 0 || /* if the feature is unset or */
                current_value == value || /* the feature is at
                                                      this value already
                                                      or */
                (current_value < 0 &&
                (-current_value != value)) /* the feature is
                                                             negatively set
                                                             to something
                                                             else */
                )
            {
                state.Set(bit_start, bit_end, value);
                return true;
            }
            return false;
        }
        throw; // for the compiler's peace of mind
    }

    typedef typename std::make_unsigned<CharType>::type UCharType;
    std::pair<UCharType, UCharType> GetFeatureValue(ccstr feature, ccstr value)const
    {
        std::pair<UCharType, UCharType> p;
        p.first = 1;
        p.second = 0;
        for (auto& f : flag_map)
        {
            if (f.first == feature)
            {
                if (value[0]) // not empty value
                {
                    p.second = 1;
                    for (auto& v : f.second)
                    {
                        if (v == value)
                        {
                            return p;
                        }
                        ++p.second;
                    }
                }
                else
                {
                    return p;
                }
            }
            ++p.first;
        }
        throw MyError("Oops!");
    }
    void CalculateOffsets()
    {
        offsets.clear();
        unsigned char bits = 0;
        if (flag_map.size() > std::numeric_limits<UCharType>::max())
        {
            throw MyError("There are ", flag_map.size(), " flag diacritic features, "
                "which is more than what a single input character can hold!");
        }
        for (const auto& flag : flag_map)
        {
            if (flag.second.size() + 1 > std::numeric_limits<UCharType>::max())
            {
                throw MyError("The diacritic flag \"", flag.first, "\" has ", flag.second.size() + 1, " possible values, "
                    "which is more than what a single input character can hold!");
            }
            const unsigned char required_flag_bits = IntLog2<unsigned char, size_t>(2 * (flag.second.size() + 1));
            offsets.emplace_back(bits);
            bits += required_flag_bits;
        }
        if (bits > sizeof(State) * CHAR_BIT)
            throw MyError("Storing the state of the flag diacritics requires ", bits,
                " bits, but FlagDiacritics::State is only ", sizeof(State) * CHAR_BIT, " bits wide!");
        offsets.emplace_back(bits);
    }
private:
    std::unordered_map<string, std::unordered_set<string>> flag_map;
    std::vector<unsigned char> offsets;
};

}
