#ifndef STAN_MATH_PRIM_SCAL_META_COMPLEX_STEP_HPP
#define STAN_MATH_PRIM_SCAL_META_COMPLEX_STEP_HPP

#include <complex>
#include <limits>
#include <cmath>

namespace stan {
namespace math {

template <class F> 
// Complex step differentiation as described at: 
// https://blogs.mathworks.com/cleve/2013/10/14/complex-step-differentiation/
double complex_step(F f, double x) {
    typedef std::complex<double> Complex;
    const double step = 
        std::sqrt(std::numeric_limits<double>::epsilon());
    Complex x_step(x, step);
    Complex f_x_step = f(x_step);
    return f_x_step.imag() / step;
}

}  // namespace math
}  // namespace stan

#endif
