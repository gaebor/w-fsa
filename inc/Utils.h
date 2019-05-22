#pragma once

#include <stdarg.h>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <string>
#include <unordered_map>
#include <map>
#include <ostream>
#include <istream>
#include <sstream>
#include <vector>
#include <utility>
#include <list>
#include <algorithm>
#include <functional>
#include <exception>
#include <thread>

#include "Isatty.h"

#include "mkl_types.h"

std::pair<const char*, char> GetWord(char*& input, const char* separator= " ");

void PrintFixedWidth(FILE* out, double x, int width = 7);

typedef const char* CStr;

struct StrEq
{
    bool operator()(const CStr& a, const CStr& b)const;
    static const StrEq streq;
};

struct StrLess
{
    bool operator()(const CStr& a, const CStr& b)const;
};

// http://sandbox.hlt.bme.hu/~gaebor/STLdoc/VS2017/xstddef.html#a4a7a1c562fda5e9e42d80a12cec4a539
struct StrHash
{
    size_t operator()(const CStr& s) const;
};

template<class T>
class Counter : public std::unordered_map<T, size_t>
{
public:
    size_t operator[](const T& key)
    {
        return this->emplace(key, this->size()).first->second;
    }
};

template<class T>
struct Keyed : std::unordered_map<CStr, T, StrHash, StrEq>
{
};

bool IsEmpty(const char* s);

std::string ToStr();

template<typename Arg1, typename... Args>
std::string ToStr(const Arg1& arg1, const Args&... args)
{
    std::ostringstream oss;
    oss << arg1;
    return oss.str() + ToStr(args...);
}
template<typename Arg1>
std::string ToStr(const Arg1& arg1)
{
    std::ostringstream oss;
    oss << arg1;
    return oss.str();
}

class MyError : public std::exception
{
public:
    template<typename ...Args>
    explicit MyError(const Args&... args)
        : msg(ToStr(args...))
    {
    }
    virtual const char* what() const throw();
private:
    std::string msg;
};

bool ReadContent(FILE* input, std::vector<char>& content);

void PrintCooMtx(std::ostream& os,
    const std::vector<double>& data,
    const std::vector<MKL_INT>& rows,
    const std::vector<MKL_INT>& cols);

void PrintCsrMtx(std::ostream& os,
    const std::vector<double>& data,
    const std::vector<MKL_INT>& rows,
    const std::vector<MKL_INT>& cols,
    const std::vector<double>& rhs = std::vector<double>());

void ReadCsrMtx(std::istream& is,
    std::vector<double>& data,
    std::vector<MKL_INT>& rows,
    std::vector<MKL_INT>& cols);

void WriteCsrMtx(std::ostream& os,
    const std::vector<double>& data,
    const std::vector<MKL_INT>& rows,
    const std::vector<MKL_INT>& cols);

bool ContainsPrefix(const CStr& word, const CStr& prefix);

double LogFactorial(size_t d);
double LogSimplexVolume(size_t d);

template<typename T, typename F, typename Func, typename ...Args>
void ProgressIndicator(T begin, T* p, F factor, const char* fmt,
                        const Func& f, Args&&... args)
{
    bool run = true;
    std::thread progress;
    if (IsStderrTty())
    {
        progress = std::thread([&]()
        {
            bool written = false;
            while (run)
            {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                if (run)
                {   
                    written = true;
                    fprintf(stderr, fmt, size_t(*p - begin)*factor);
                    fflush(stderr);
                }
            }
            if (written)
            {
                fprintf(stderr, fmt, size_t(*p - begin)*factor);
                fflush(stderr);
            }
        });
    }
    try
    {
        f(args...);
    }
    catch (const MyError& e)
    {
        run = false;
        if (progress.joinable())
            progress.join();
        throw e;
    }
    run = false;
    if (progress.joinable())
        progress.join();
}

//! kind-of a vector implementation of set
template<typename Ty, typename Allocator>
typename std::vector<Ty, Allocator>::iterator SortedInsert(std::vector<Ty, Allocator>& vec, const Ty& value)
{
    auto where = std::upper_bound(vec.begin(), vec.end(), value);
    if (where == vec.end() || *where != value)
        return vec.insert(where, value);
    else
        return where;
}

template<typename Ty1, typename Ty2>
struct Lesser
{
    bool operator()(const std::pair<Ty1, Ty2>& one, const std::pair<Ty1, Ty2>& other)const
    {
        return one.first < other.first;
    }
};

//! kind-of a vector implementation of map
template<typename Key, typename Value, typename Allocator>
Value& SortedInsert(std::vector<std::pair<Key, Value>, Allocator>& vec, const Key& key)
{
    static std::pair<Key, Value> pair;
    pair.first = key;
    auto where = std::upper_bound(vec.begin(), vec.end(), pair, Lesser<Key, Value>());
    if (where == vec.end() || where->first != key)
        return vec.emplace(where, pair)->second;
    else
        return where->second;
}

template<typename Value>
Value& GetCoord(const std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols, std::vector<Value>& data, MKL_INT i, MKL_INT j)
{
    auto where = std::lower_bound(cols.data() + rows[i], cols.data() + rows[i + 1], j);
    return data[where - cols.data()];
}
