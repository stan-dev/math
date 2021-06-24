#ifndef STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP
#define STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {


/**
 * Return the real part of the complex argument.
 *
 * @return a complex type with real and imaginary set to 0
 */
inline std::complex<double> to_complex(){
    return std::complex<double>(0, 0);
} 


/**
 * Return the real part of the complex argument.
 *
 * @tparam T value type of argument
 * @param[in] re real value argument
 * @param[in] im imaginary value argument
 * @return complex type with real value set to re and imaginary
 *         set im or 0 if no second argument is supplied
 */
template <typename T>
inline std::complex<T> to_complex(const T& re, const T& im = 0){
    return std::complex<T>(re, im);
} 



}  // namespace math
}  // namespace stan

#endif
