#ifndef STAN_MATH_PRIM_FUN_IMAG_HPP
#define STAN_MATH_PRIM_FUN_IMAG_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the imaginary part of the complex argument.
 *
 * @tparam T value type of argument
 * @param[in] z argument
 * @return imaginary part of argument
 */
template <typename T, require_autodiff_t<T>>
T imag(const std::complex<T>& z) {
  return z.imag();
}

}  // namespace math
}  // namespace stan

#endif
