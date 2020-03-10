#include "Utils.h"

#include <algorithm>
#include <vector>

// #include "mkl_dss.h"
#define _USE_MATH_DEFINES 
#include <math.h>

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
    // GO TO HELL
    const int magnitude = (x == 0) ? 0 : (int)std::floor(std::log10(std::abs(x)));
    if (magnitude <= width - 2 && magnitude >= 0)
    {
        if (std::floor(x) == x)
            fprintf(out, "%*.0f", width, x);
        else
            fprintf(out, "%*.*f", width, std::max(0, width - 3 - magnitude), x);
    }
    else if (-4 < magnitude && magnitude < 0)
        fprintf(out, "%*.*f", width, width - 3, x);
    else
        fprintf(out, "%*.*e", width, width - 7, x);
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

//void PrintCooMtx(std::ostream& os,
//    const std::vector<double>& data,
//    const std::vector<MKL_INT>& rows,
//    const std::vector<MKL_INT>& cols)
//{
//    os.setf(std::ios::fixed);
//
//    const auto max_row = (*std::max_element(rows.begin(), rows.end())) + 1;
//    const auto oldprec = os.precision(4);
//    for (MKL_INT r = 0; r < max_row; ++r)
//    {
//        for (MKL_INT i = cols.size()-1; i >= 0; --i)
//            if (rows[i] == r)
//            {
//                os << '\r';
//                for (MKL_INT j = 0; j < cols[i]; ++j) os << "       ";
//                os << data[i];
//            }
//        os << '\n';
//    }
//    os.precision(oldprec);
//    os.unsetf(std::ios::fixed);
//}
//
//void PrintCsrMtx(FILE* f, const std::vector<double>& data,
//    const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols,
//    const std::vector<double>& rhs)
//{
//    const MKL_INT width = rows.size() - 1;
//    MKL_INT prev;
//    for (size_t i = 0; i < rows.size()-1; ++i)
//    {
//        prev = -1;
//        for (MKL_INT j = rows[i]; j < rows[i + 1]; ++j)
//        {
//            for (; prev < cols[j]-1; ++prev) fputs("        ", f);
//            
//            prev = cols[j];
//            PrintFixedWidth(f, data[j], 7);
//            fputs(" ", f);
//        }
//        if (rhs.size() > i)
//        {
//            for (; prev < width -1 ; ++prev) fputs("        ", f);
//            fputs("|", f);
//            PrintFixedWidth(f, rhs[i], 7);
//        }
//        fputs("\n", f);
//    }
//}
//
//void ReadCsrMtx(std::istream & is, std::vector<double>& data, std::vector<MKL_INT>& rows, std::vector<MKL_INT>& cols)
//{
//    data.clear(); rows.clear(); cols.clear();
//    std::string line;
//    while (std::getline(is, line))
//    {
//        rows.push_back(cols.size());
//        std::istringstream iss(line);
//        do
//        {
//            cols.emplace_back();
//            data.emplace_back();
//            iss >> cols.back() >> data.back();
//        } while (iss);
//        cols.pop_back();
//        data.pop_back();
//    }
//    rows.push_back(cols.size());
//}
//
//void WriteCsrMtx(std::ostream & os, const std::vector<double>& data, const std::vector<MKL_INT>& rows, const std::vector<MKL_INT>& cols)
//{
//    for (size_t rowi = 0; rowi < rows.size() - 1; ++rowi)
//    {
//        for (MKL_INT colj = rows[rowi]; colj < rows[rowi + 1]; ++colj)
//        {
//            os << cols[colj] << ' ' << data[colj] << ' ';
//        }
//        os << '\n';
//    }
//}

bool ContainsPrefix(const CStr& word, const CStr& prefix)
{
    return strncmp(word, prefix, strlen(prefix)) == 0;
}

const char* ContainsPrefix2(const char* word, const char* prefix)
{
    while (*word == *prefix && *word != '\0')
    {
        ++word;
        ++prefix;
    }
    return (*prefix == '\0') ? word : nullptr;
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

double mxlogx(double x)
{
    if (x > 0)
        return x * (-std::log(x));
    else if (x == 0)
        return 0;
    else
        return std::numeric_limits<decltype(x)>::infinity();
}

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

const char* MyError::what() const throw()
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

//double RealSymmetricLogDet(const MKL_INT* const Hrow, const MKL_INT n,
//                       const MKL_INT* const Hcol, const MKL_INT nnz,
//                       const double* const Hdata, bool reorder)
//{
//    if (n == nnz)
//    {
//        double result = 0;
//        for (MKL_INT i = 0; i < n; ++i)
//        {
//            if (Hdata[i] > 0)
//                result += std::log(Hdata[i]);
//            else
//                return atof("inf");
//        }
//        return result;
//    }
//    std::vector<MKL_INT> perm;
//    MKL_INT solver_opt = MKL_DSS_ZERO_BASED_INDEXING +
//        MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR +
//        MKL_DSS_REFINEMENT_ON;
//    
//    MKL_INT error;
//    DssSolverHandler solver(solver_opt);
//    auto solver_handle = solver.GetHandler();
//    
//    solver_opt = MKL_DSS_SYMMETRIC;
//    if ((error = dss_define_structure(solver_handle, solver_opt, Hrow, n, n, Hcol, nnz)) != MKL_DSS_SUCCESS)
//    {
//        throw MyError("Unable to define structure for RealSymmetricLogDet! Error code: ", error);
//    }
//    solver_opt = reorder ? MKL_DSS_GET_ORDER : MKL_DSS_MY_ORDER;
//    perm.resize(n);
//    if (!reorder)
//    {   // stick to original order
//        for (MKL_INT i = 0; i < n; ++i)
//            perm[i] = i;
//    }
//    if ((error = dss_reorder(solver_handle, solver_opt, perm.data())) != MKL_DSS_SUCCESS)
//    {
//        throw MyError("Unable to find reorder for RealSymmetricLogDet! Error code: ", error);
//    }
//    solver_opt = MKL_DSS_INDEFINITE;
//    if ((error = dss_factor_real(solver_handle, solver_opt, Hdata)) != MKL_DSS_SUCCESS)
//    {
//        throw MyError("Unable to factor for RealSymmetricLogDet! Error code: ", error);
//    }
//    solver_opt = MKL_DSS_DEFAULTS;
//    double det[2];
//    if ((error = dss_statistics(solver_handle, solver_opt, "Determinant", det)) != MKL_DSS_SUCCESS)
//    {
//        throw MyError("Unable to get determinant for RealSymmetricLogDet! Error code: ", error);
//    }
//    else
//    {
//        const double& det_pow = det[0], &det_base = det[1];
//        if (det_base <= 0.0)
//            return atof("inf");
//        return std::log(det_base) + det_pow * M_LN10;
//    }
//}
//
//DssSolverHandler::DssSolverHandler(MKL_INT solver_opt)
//    : solver_handler(NULL)
//{
//    MKL_INT error;
//    if ((error = dss_create(solver_handler, solver_opt)) != MKL_DSS_SUCCESS)
//    {
//        throw MyError("Failed to create DSS! Error code: ", error);
//    }
//}
//
//DssSolverHandler::~DssSolverHandler()
//{
//    if (solver_handler)
//    {
//        MKL_INT solver_opt = MKL_DSS_MSG_LVL_WARNING + MKL_DSS_TERM_LVL_ERROR;
//        dss_delete(solver_handler, solver_opt);
//        //if (error != MKL_DSS_SUCCESS)
//        //    throw MyError("Cannot delete solver! Error: ", error);
//        solver_handler = 0;
//    }
//}
//
//_MKL_DSS_HANDLE_t DssSolverHandler::GetHandler() const
//{
//    return solver_handler;
//}
//
//SparseMtxHandle::SparseMtxHandle()
//{
//    desc.type = SPARSE_MATRIX_TYPE_GENERAL;
//}
//
//SparseMtxHandle::SparseMtxHandle(MKL_INT n, MKL_INT k, MKL_INT * rows, MKL_INT * cols, double * x)
//{
//    desc.type = SPARSE_MATRIX_TYPE_GENERAL;
//    Init(n, k, rows, cols, x);
//}
//
//void SparseMtxHandle::Init(MKL_INT n, MKL_INT k, MKL_INT * rows, MKL_INT * cols, double * x)
//{
//    auto status = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, n, k, rows, rows + 1, cols, x);
//    if (SPARSE_STATUS_SUCCESS != status)
//    {
//        A = nullptr;
//        throw MyError("mkl_sparse_d_create_csr failed with ", status);
//    }
//}
//
//SparseMtxHandle::~SparseMtxHandle()
//{
//    if (A)
//        mkl_sparse_destroy(A);
//}
//
//void SparseMtxHandle::dot(const double * x, double * y, sparse_operation_t op, double alpha, double beta) const
//{
//    auto status = mkl_sparse_d_mv(op, alpha, A, desc, x, beta, y);
//    if (status != SPARSE_STATUS_SUCCESS)
//    {
//        throw MyError("mkl_sparse_d_mv failed with ", status);
//    }
//}
