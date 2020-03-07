#include "lookup.h"

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>

Context* self;

static inline void push_back(unsigned int i)
{
    if (self->size == self->allocated)
    {
        self->path = realloc(self->path, ((self->allocated) *= 2)*sizeof(unsigned int));
    }
    (self->path)[(self->size)++] = i;
}

static inline void pop_back()
{
    if (self->size > 0)
        --(self->size);
}

static inline int str_ends(unsigned int x)
{
    for (size_t i = 0; i < sizeof(unsigned int); ++i)
    {
        if (((const char*)(&x))[i] == '\0')
            return 1;
    }
    return 0;
}

static inline const char* GetNextUtf8Character(const char* word)
{
    if ((128 & *word) == '\0')
    {
        return word + 1;
    }
    else if ((32 & *word) == '\0')
    {
        return word + 2;
    }
    else if ((16 & *word) == '\0')
    {
        return word + 3;
    }
    else
    {
        return word + 4;
    }
}

static inline const char* ContainsPrefix2(const char* word, const char* prefix)
{
    while (*word == *prefix && *word != '\0')
    {
        ++word;
        ++prefix;
    }
    return (*prefix == '\0') ? word : NULL;
}

void lookup(const char* s, unsigned int beg, const unsigned int end)
{
    size_t i;
    const char* next;

    for (; beg < end; ++beg)
    {   // try outgoing edges
        const unsigned int id = self->table[beg];
        const unsigned int to_beg = self->table[beg + 1];
        const unsigned int to_end = self->table[beg + 2];
        const char* input = (const char*)(self->table + beg + 3);

        if (to_end == UINT_MAX || to_beg == UINT_MAX)
        {   //final state
            if (*s == '\0')
            {   // that's a result
                push_back(id);
                for (i = 0; i < self->size; ++i)
                {
                    fprintf(stdout, "%u ", (self->path)[i]);
                }
                putc('\n', stdout);
                self->has_analysis = 1;
                pop_back();
            }
        }
        else if (strcmp(input, "@_UNKNOWN_SYMBOL_@") == 0 || strcmp(input, "@_IDENTITY_SYMBOL_@") == 0)
        {   // consume one character
            // TODO this can be hastened is the special symbols are shorter!
            push_back(id);
            lookup(GetNextUtf8Character(s), to_beg, to_end);
        }
        // 
        // epsilon is handled with a simple empty string
        // 
        else if (input[0] == '@' && input[2] == '.' &&
            (input[1] == 'P' || input[1] == 'N' || input[1] == 'D' || input[1] == 'R' || input[1] == 'C' || input[1] == 'U') &&
            input[3] != '\0')
        {
            // FlagDiacritics::State flagsate = fd_state;
            // if (FlagDiacritics::Apply(input, fd_state))
            {
                push_back(id);
                lookup(s, to_beg, to_end);
            }
            // fd_state = flagsate;
        }
        else if ((next = ContainsPrefix2(s, input)) != NULL)
        {   // a lead to follow
            push_back(id);
            lookup(next, to_beg, to_end);
        }

        beg += 3;
        while (str_ends((self->table)[beg]) == 0)
        {
            ++beg;
        }
    }

    pop_back();
}
