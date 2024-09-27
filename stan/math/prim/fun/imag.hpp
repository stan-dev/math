#ifndef STAN_MATH_PRIM_FUN_IMAG_HPP
#define STAN_MATH_PRIM_FUN_IMAG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/complex_base.hpp>

namespace stan {
namespace math {

/**
 * Return the imaginary component of the complex argument.
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose imaginary component is extracted
 * @return imaginary component of argument
 */
template <typename T, require_autodiff_t<T>>
T imag(const stan::math::complex<T>& z) {
  return z.imag();
}

}  // namespace math
}  // namespace stan

#endif
