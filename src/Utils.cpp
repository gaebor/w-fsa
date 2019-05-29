#include "Utils.h"

#include <algorithm>
#include <vector>

bool matches(const char * str, const char * pattern1)
{
    return (strcmp(pattern1, str) == 0);
}

bool matches(const char *)
{
    return false;
}

std::pair<const char*, char> GetWord(char*& input, const char* separator)
{
    const char* result = input;
    const char* sep = separator;
    enum {
        IN_SEP,
        IN_WORD,
    } s = IN_WORD;

    while (*input)
    {
        switch (s)
        {
        case IN_WORD:
            if (*input == *sep)
            {
                s = IN_SEP;
                ++sep;
            }
            else if (*input == '\n')
            {
                *(input++) = '\0';
                return { result, '\n' };
            }
            break;
        case IN_SEP:
            if (*input == *sep)
            {
                // consume separator
                ++sep;
            }
            else if (*input == '\n' && *sep == '\0')
            {   // separator successfully consumed but also end of line appeared
                for (char* w = input - (sep - separator); w < input; ++w)
                    *w = '\0';
                return { result, '\n' };
            }
            else if (*input == '\n')
            {   // end of line appeared inside a separator
                *(input++) = '\0';
                return { result, '\n' };
            }
            else if (*sep == '\0')
            {
                // separator successfully consumed inside a line
                for (char* w = input - (sep - separator); w < input; ++w)
                    *w = '\0';
                return { result, *(sep - 1) };
            }
            else
            {
                // not a separator after all, reset
                sep = separator;
                s = IN_WORD;
            }
            break;
        }
        ++input;
    }
    return { result, '\0'};
}

void PrintFixedWidth(FILE * out, double x, int width)
{
    fprintf(out, "%*.*e", width, std::max(0, width - 7), x);
    //static char fmt[10];
    //snprintf(fmt, 9, "%%-%d.*g", width);
    //if (x == 0)
    //    fprintf(out, "%g", x);
    //else
    //{
    //    const double d = floor(abs(log10(abs(x))));
    //    if (x < 0)
    //        fprintf(out, fmt, std::max(std::min(width - 3 - (int)d, width - 3), 0), x);
    //    else
    //        fprintf(out, fmt, std::max(std::min(width - 2 - (int)d, width - 2), 0), x);
    //}
}

bool IsEmpty(const char * s)
{
    return s[0] == '\0';
}

std::string ToStr()
{
    return std::string();
}

bool ReadContent(FILE * input, std::vector<char>& content)
{
    if (input)
    {
        auto begin = ftell(input);
        if (fseek(input, 0, SEEK_END) != 0)
            return false;
        auto length = ftell(input);
        if (length == -1L)
            return false;
        length -= begin;
        if (fseek(input, begin, SEEK_SET) != 0)
            return false;

        content.resize(length);
        content.resize(fread(content.data(), 1, length, input));
        content.emplace_back('\0');
    }
    else
        return false;

    return true;
}

void PrintCooMtx(std::ostream & os,
    const std::vector<double>& data,
    const std::vector<MKL_INT>& rows,
    const std::vector<MKL_INT>& cols)
{
    os.setf(std::ios::fixed);

    const auto max_row = (*std::max_element(rows.begin(), rows.end())) + 1;
    const auto oldprec = os.precision(4);
    for (MKL_INT r = 0; r < max_row; ++r)
    {
        for (MKL_INT i = cols.size()-1; i >= 0; --i)
            if (rows[i] == r)
            {
                os << '\r';
                for (MKL_INT j = 0; j < cols[i]; ++j) os << "       ";
                os << data[i];
            }
        os << '\n';
    }
    os.precision(oldprec);
    os.unsetf(std::ios::fixed);
}

void PrintCsrMtx(std::ostream & os, const std::vector<double>& data,
    const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols,
    const std::vector<double>& rhs)
{
    const MKL_INT width = rows.size() - 1;
    MKL_INT prev;
    os.setf(std::ios::fixed);
    const auto oldprec = os.precision(4);
    for (size_t i = 0; i < rows.size()-1; ++i)
    {
        prev = -1;
        for (MKL_INT j = rows[i]; j < rows[i + 1]; ++j)
        {
            for (; prev < cols[j]-1; ++prev) os << "       ";
            
            prev = cols[j];
            os << data[j] << ' ';
        }
        if (rhs.size() > i)
        {
            for (; prev < width -1 ; ++prev) os << "       ";
            os << '|' << rhs[i];
        }
        os << "\n";
    }
    os.precision(oldprec);
    os.unsetf(std::ios::fixed);
}

void ReadCsrMtx(std::istream & is, std::vector<double>& data, std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols)
{
    data.clear(); rows.clear(); cols.clear();
    std::string line;
    while (std::getline(is, line))
    {
        rows.push_back(cols.size());
        std::istringstream iss(line);
        do
        {
            cols.emplace_back();
            data.emplace_back();
            iss >> cols.back() >> data.back();
        } while (iss);
        cols.pop_back();
        data.pop_back();
    }
    rows.push_back(cols.size());
}

void WriteCsrMtx(std::ostream & os, const std::vector<double>& data, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols)
{
    for (size_t rowi = 0; rowi < rows.size() - 1; ++rowi)
    {
        for (MKL_INT colj = rows[rowi]; colj < rows[rowi + 1]; ++colj)
        {
            os << cols[colj] << ' ' << data[colj] << ' ';
        }
        os << '\n';
    }
}

bool ContainsPrefix(const CStr& word, const CStr& prefix)
{
    return strncmp(word, prefix, strlen(prefix)) == 0;
}

double LogSimplexVolume(size_t d)
{
    if (d > 0)
        return 0.5*std::log(d) - LogFactorial(d - 1);
    else
        return 0;
}

struct LogFactorialMemory : std::vector<double>
{
    LogFactorialMemory() : std::vector<double>(2, 0) {}
};

double LogFactorial(size_t d)
{
    static LogFactorialMemory memory;
    
    if (memory.size() > d)
        return memory[d];
    else
    {
        double result = memory.back();
        memory.reserve(d + 1);
        for (auto i = memory.size(); i <= d; ++i)
        {
            result += std::log(i);
            memory.emplace_back(result);
        }
        return result;
    }
}

const char * MyError::what() const throw()
{
    return msg.c_str();
}

bool StrEq::operator()(const CStr& a, const CStr& b) const
{
    return strcmp(a, b) == 0;
}
bool StrLess::operator()(const CStr& a, const CStr& b) const
{
    return strcmp(a, b) < 0;
}

const StrEq StrEq::streq;

struct FNV
{
    static const size_t offset = ((sizeof(size_t) == 8) ? 14695981039346656037ULL : 2166136261U);
    static const size_t prime = ((sizeof(size_t) == 8) ? 1099511628211ULL : 16777619U);
};

size_t StrHash::operator()(const CStr& s) const
{
    size_t result(FNV::offset);
    for (auto p = s; *p; ++p)
    {
        result ^= (size_t)(*p);
        result *= FNV::prime;
    }
    return result;
}
