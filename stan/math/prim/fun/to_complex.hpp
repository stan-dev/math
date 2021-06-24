#ifndef STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP
#define STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {



//Create a complex type from a real part and an imaginary part
template <typename T>
std::complex<T> to_complex(const T& re, const T& im = 0){
    return std::complex<T>(re, im);
} 


}  // namespace math
}  // namespace stan

#endif
