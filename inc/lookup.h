#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
    unsigned int* table;
    unsigned int* path;
    size_t size, allocated;
    char has_analysis;
} Context;

extern Context* self;

void lookup(const char* s, unsigned int beg, const unsigned int end);

#ifdef __cplusplus
}
#endif
