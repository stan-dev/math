#ifndef STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP
#define STAN_MATH_PRIM_FUN_TO_COMPLEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return a complex value from a real component and an imaginary component.
 * Default values for both components is 0.
 *
 * @tparam T type of real component
 * @tparam S type of imaginary component
 * @param[in] re real component (default = 0)
 * @param[in] im imaginary component (default = 0)
 * @return complex value with specified real and imaginary components
 */
template <typename T = double, typename S = double>
constexpr inline std::complex<stan::real_return_t<T, S>> to_complex(
    const T& re = 0, const S& im = 0) {
  return std::complex<stan::real_return_t<T, S>>(re, im);
}

}  // namespace math
}  // namespace stan

#endif
