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
template<typename Value, typename Allocator>
Value& SortedInsert(std::vector<Value, Allocator>& vec, const Value& value)
{
    auto where = std::lower_bound(vec.begin(), vec.end(), value);
    if (where == vec.end() || *where != value)
        return *vec.emplace(where, value);
    else
        return *where;
}

template<typename Ty1, typename Ty2>
struct Lesser
{
    bool operator()(const std::pair<Ty1, Ty2>& one, const std::pair<Ty1, Ty2>& other)const
    {
        return one.first < other.first;
    }
};

template<class Vec1, class Iterator2, class Compare>
void Intersect(Vec1& v1, Iterator2 first2, Iterator2 last2, const Compare& comp)
{
    auto first1 = v1.begin();
    auto last1 = v1.end();
    
    while (first1 != last1 && first2 != last2)
    {
        if (comp.less(*first1, first2))
        {   // remove from first
            first1 = v1.erase(first1);
            // iterator validity!
            last1 = v1.end();
        }
        else if (comp.greater(*first1, first2))
        {   // skip from second
            ++first2;
        }
        else
        {   // leave it
            ++first1;
            ++first2;
        }
    }
    // remove the rest, if any
    v1.erase(first1, last1);
}
template<class Vec1, class Iterator2, class Compare>
void Subtract(Vec1& v1, Iterator2 first2, Iterator2 last2, const Compare& comp)
{
    auto first1 = v1.begin();
    auto last1 = v1.end();

    while (first1 != last1 && first2 != last2)
    {
        if (comp(*first1, *first2))
        {   // keep it
            ++first1;
        }
        else if (comp(*first2, *first1))
        {   // skip from second
            ++first2;
        }
        else
        {   // remove from first
            first1 = v1.erase(first1);
            // iterator validity!
            last1 = v1.end();
            ++first2;
        }
    }
}

template<class Vec1, class Iterator2>
void Subtract(Vec1& v1, Iterator2 first2, Iterator2 last2)
{
    Subtract<Vec1, Iterator2, std::less<typename Vec1::value_type>>(v1, first2, last2, std::less<typename Vec1::value_type>());
}

template<class Vec1, class Iterator2, class Compare>
void Union(Vec1& v1, Iterator2 first2, Iterator2 last2, Compare comp)
{
    static_assert(std::is_same<typename Vec1::value_type, typename Iterator2::value_type>::value, "the two iterables should have the same value_type!");
    auto first1 = v1.begin();
    auto last1 = v1.end();

    while (first1 != last1 && first2 != last2)
    {
        if (comp(*first1, *first2))
        {
            ++first1;
        }
        else if (comp(*first2, *first1))
        {
            first1 = v1.insert(first1, *first2);
            // iterator validity!
            last1 = v1.end();
            ++first1;
            ++first2;
        }
        else
        {
            ++first1;
            ++first2;
        }
    }
    // add the rest, if any
    v1.insert(v1.end(), first2, last2);
}

template<class Vec1, class Iterator2>
void Union(Vec1& v1, Iterator2 first2, Iterator2 last2)
{
    Union<Vec1, Iterator2, std::less<typename Vec1::value_type>>(v1, first2, last2, std::less<typename Vec1::value_type>());
}

double mxlogx(double x);

//! kind-of a vector implementation of map
template<typename Key, typename Value, typename Allocator>
Value& SortedInsert(std::vector<std::pair<Key, Value>, Allocator>& vec, const Key& key)
{
    std::pair<Key, Value> pair;
    pair.first = key;
    auto where = std::lower_bound(vec.begin(), vec.end(), pair, Lesser<Key, Value>());
    if (where == vec.end() || where->first != key)
        return vec.emplace(where, pair)->second;
    else
        return where->second;
}

template<typename Value>
Value& GetCoord(const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, std::vector<Value>& data, MKL_INT i, MKL_INT j)
{
    auto where = std::lower_bound(cols.data() + rows[i], cols.data() + rows[i + 1], j);
    return data[where - cols.data()];
}

template<typename Value>
Value GetCoord2(const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols, const std::vector<Value>& data, MKL_INT i, MKL_INT j)
{
    const auto begin = cols.data() + rows[i];
    const auto end = cols.data() + rows[i + 1];
    auto where = std::lower_bound(begin, end, j);
    return (where == end || *where != j) ? 0.0 : data[where - cols.data()];
}

struct DssSolverHandler
{
    DssSolverHandler(MKL_INT solver_opt);
    ~DssSolverHandler();
    void* GetHandler()const;
private:
    void* solver_handler;
};

double RealSymmetricLogDet(const MKL_INT* const Hrow, const MKL_INT n,
                       const MKL_INT* const Hcol, const MKL_INT nnz,
                       const double* const Hdata, bool reorder=false);

//! iterator for iterating only on first value of pairs
template<class Iterator>
struct FirstIterator
{
    typedef typename Iterator::value_type::first_type value_type;
    FirstIterator(const Iterator& it) : _it(it) {}
    FirstIterator& operator++()
    {
        ++_it;
        return *this;
    }
    value_type& operator*()const
    {
        return _it->first;
    }
    bool operator==(const FirstIterator& other)const
    {
        return _it == other._it;
    }
    bool operator!=(const FirstIterator& other)const
    {
        return _it != other._it;
    }
    Iterator _it;
};

template<class Iterator>
FirstIterator<Iterator> make_first_iterator(const Iterator& it)
{
    return FirstIterator<Iterator>(it);
}

//template<class Iterable>
//void Dok2Csr(const Iterable& indices, std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols, bool include_diag = false)
//{
//    rows.clear(); cols.clear();
//    MKL_INT row = -1;
//    for (const auto& ij : indices)
//    {
//        if (ij.first > row)
//        {
//            for (MKL_IT r = row + 1; r < ij.first; ++r)
//            {   // fill in the missing diagonals
//
//            }
//            rows.emplace_back(cols.size());
//            row = ij.first;
//        }
//        Hcol.emplace_back(ij.second);
//    }
//    rows.emplace_back(cols.size());
//
//}
