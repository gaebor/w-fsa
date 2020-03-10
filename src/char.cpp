#include "atto/char.h"

namespace atto {

template<>
const typename CodeUnit<Encoding::UTF8>::type*& StepNextCharacter<Encoding::UTF8>(const typename CodeUnit<Encoding::UTF8>::type*& word)
{
    // step over continuation bytes
    do
    {
        ++word;
    }while(((*word) & 0b11000000) == 0b10000000);
    return word;
}

template<>
const typename CodeUnit<Encoding::UTF16>::type*& StepNextCharacter<Encoding::UTF16>(const typename CodeUnit<Encoding::UTF16>::type*& word)
{
    if (((*word) & 0b1111110000000000) == 0b1101100000000000)
        word += 2;
    else
        ++word;
    return word;
}

}
