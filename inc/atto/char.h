#pragma once

#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace atto {

enum Encoding
{
    ASCII,
    CP,
    OCTET = CP,
    UTF8,
    UCS2,
    UTF16,
    UTF32
};

template<Encoding e>
struct CodeUnit
{
    typedef char type;
    static const size_t size = sizeof(type);
    type c;
};

template<>
struct CodeUnit<Encoding::UCS2>
{
#ifdef WIN32
    typedef wchar_t type;
#else
    typedef std::uint16_t type;
#endif
    static const size_t size = sizeof(type);
    type c;
};

// see https://en.cppreference.com/w/cpp/language/types#Character_types
template<>
struct CodeUnit<Encoding::UTF16>
{
#ifdef WIN32
    typedef wchar_t type;
#else
    typedef std::uint16_t type;
#endif
    static const size_t size = sizeof(type);
    type c;
};

// see https://en.cppreference.com/w/cpp/language/types#Character_types
template<>
struct CodeUnit<Encoding::UTF32>
{
#ifdef WIN32
    typedef std::uint32_t type;
#else
    typedef wchar_t type;
#endif
    static const size_t size = sizeof(type);
    type c;
};

template<Encoding e>
const typename CodeUnit<e>::type*& StepNextCharacter(const typename CodeUnit<e>::type*& word)
{
    return ++word;
}

template<>
const typename CodeUnit<Encoding::UTF8>::type*& StepNextCharacter<Encoding::UTF8>(const typename CodeUnit<Encoding::UTF8>::type*& word);

template<>
const typename CodeUnit<Encoding::UTF16>::type*& StepNextCharacter<Encoding::UTF16>(const typename CodeUnit<Encoding::UTF16>::type*& word);

template<Encoding e>
const typename CodeUnit<e>::type* GetNextCharacter(const typename CodeUnit<e>::type* word)
{
    return StepNextCharacter<e>(word);
}

template<class CharType, class StorageType>
bool StrEnds(StorageType word)
{
    static_assert(std::is_integral<CharType>::value, "");
    for (size_t i = 0; i < sizeof(StorageType); ++i)
    {
        if (((const CharType*)(&word))[i] == (CharType)0)
            return true;
    }
    return false;
}

template<Encoding e, class StorageType>
bool StrEnds(StorageType word)
{
    return StrEnds<typename CodeUnit<e>::type>(word);
}

template<class CharType>
bool StrEqual(const CharType* word1, const CharType* word2)
{
    static_assert(std::is_integral<CharType>::value, "");
    while (*word1 == *word2 && *word1 != (CharType)0)
    {
        ++word1;
        ++word2;
    }
    return *word1 == (CharType)0 && *word2 == (CharType)0;
}

template<Encoding e>
bool StrEqual(const typename CodeUnit<e>::type* word1, const typename CodeUnit<e>::type* word2)
{
    return StrEqual<typename CodeUnit<e>::type>(word1, word2);
}

template<class CharType>
const CharType* StrPrefix(const CharType* word, const CharType* prefix)
{
    static_assert(std::is_integral<CharType>::value, "");
    while (*word == *prefix && *word != (CharType)0)
    {
        ++word;
        ++prefix;
    }
    return (*prefix == (CharType)0) ? word : nullptr;
}

template<Encoding e>
const typename CodeUnit<e>::type* StrPrefix(const typename CodeUnit<e>::type* word, const typename CodeUnit<e>::type* prefix)
{
    return StrPrefix<typename CodeUnit<e>::type>(word, prefix);
}

}
