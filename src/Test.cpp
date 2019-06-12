
#include "Test.h"

#include "mkl.h"

#include <chrono>
#include <vector>
#include <cmath>
#include <numeric>

//! simple class to measure time, not thread safe!
/*! from https://github.com/gaebor/AsyncQueue
*/
template<class ClockTy = std::chrono::steady_clock>
class Clock
{
public:
    //! marks instantiation time
    Clock()
    {
        Tick();
    }
    //! marks current time
    void Tick()
    {
        _timePoint = MyClock::now();
    }
    //! returns time between former time mark and now
    /*!
    @return time since last Tick, construction.
    */
    double Tack()
    {
        auto const now = MyClock::now();
        const auto elapsed = _frequency * (now - _timePoint).count();
        return elapsed;
    }
    ~Clock() {}
private:
    typedef ClockTy MyClock;
    typedef typename MyClock::duration MyDuration;
    typename MyClock::time_point _timePoint;
    static constexpr double _frequency = (double)MyDuration::period::num / MyDuration::period::den;
};

void test(std::ostream& os, size_t epoch, size_t nn, size_t tests)
{
    Clock<> clock;
    std::vector<double> out, x;

    os << "n" <<
        ((tests & 1) ? "\texp(x)" : "") <<
        ((tests & 2) ? "\tmkl" : "") <<
        ((tests & 4) ? "\tsum(x)" : "") <<
        ((tests & 8) ? "\tx.1" : "") << 
        ((tests & 16) ? "\tx=0" : "") << 
        ((tests & 32) ? "\tx*=0" : "") << std::endl;

    os.precision(3);
    os << std::scientific;

    for (size_t m = 0; m < 32; ++m)
    {
        if (0 == (nn & ((size_t)1 << m)))
            continue;
        const auto n = (MKL_INT)std::pow(10, m);
        
        os << "10^" << m << '\t';
        out.resize(n);
        double sum = 0.0;

        x.assign(n, 3.14);

        if (tests & 1)
        {
            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                for (MKL_INT i = 0; i < n; ++i)
                    out[i] = std::exp(x[i]);

            os << (clock.Tack() / epoch) << '\t';
            os.flush();
        }
        if (tests & 2)
        {
            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                vdExp(n, x.data(), out.data());

            os << (clock.Tack() / epoch) << '\t';
            os.flush();
        }
        if (tests & 4)
        {
            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                sum = std::accumulate(x.begin(), x.end(), 0.0);

            os << (clock.Tack() / epoch) << '\t';
            os.flush();
        }
        if (tests & 8)
        {
            out.assign(n, 1.0);

            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                cblas_ddot(n, x.data(), 1, out.data(), 1);
            os << (clock.Tack() / epoch) << "\t";
            os.flush();
        }
        if (tests & 16)
        {
            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                std::fill(x.begin(), x.end(), 0);
            os << (clock.Tack() / epoch) << "\t";
            os.flush();
        }
        if (tests & 32)
        {
            clock.Tick();
            for (size_t e = 0; e < epoch; ++e)
                cblas_dscal(n, 0.0, x.data(), 1);
            os << (clock.Tack() / epoch) << "\t";
            os.flush();
        }
        os << "(" << sum << ")" << std::endl;
    }

}
