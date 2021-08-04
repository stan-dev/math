#ifndef STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP
#define STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <type_traits>
#include <complex>

namespace stan {
namespace math {



/**
 * Return a complex type from a real part and an imaginary part.
 * Default values for both parts is 0.
 *
 * @tparam T type of real and
 * @param re real element
 * @param im imaginary element
 * @return Complex type
 */
template <typename T = double, typename S = double>
inline std::complex<T> to_complex(const T& re = 0, const S& im = 0){
    return std::complex<T>(re, im);
}

/**
 * Return a complex type from a integer real part and an imaginary part.
 *
 * @tparam T type of real and
 * @param re real element
 * @param im imaginary element
 * @return Complex type
 *
 * Used for type casting from int to double for real component.
 * Otherwise arrays of complex don't work with ints.
 */
template <typename S = double>
inline std::complex<double> to_complex(const int& re, const S& im = 0){
    return std::complex<double>(re, im);
}

}  // namespace math
}  // namespace stan

#endif
