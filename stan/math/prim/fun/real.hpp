#ifndef STAN_MATH_PRIM_FUN_REAL_HPP
#define STAN_MATH_PRIM_FUN_REAL_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the real component of the complex argument.
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose real component is extracted
 * @return real component of argument
 */
template <typename T, require_autodiff_t<T>>
inline T real(const std::complex<T>& z) {
  return z.real();
}

/**
 * Return the real component of the complex argument.
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose real component is extracted
 * @return real component of argument
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline T& real(std::complex<T>& z) {
  return reinterpret_cast<T(&)[2]>(z)[0];
}

/**
 * Return the real component of the complex argument.
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose real component is extracted
 * @return real component of argument
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline T& real(std::complex<T>&& z) {
  return reinterpret_cast<T(&)[2]>(z)[0];
}


}  // namespace math
}  // namespace stan

#endif
