#ifndef STAN_MATH_PRIM_FUN_IMAG_HPP
#define STAN_MATH_PRIM_FUN_IMAG_HPP

#include <stan/math/prim/meta.hpp>
#include <complex>

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
inline T imag(const std::complex<T>& z) {
  return z.imag();
}

/**
 * Return the imaginary component of the complex argument.
 * For any object z of type std::complex<T>, reinterpret_cast<T(&)[2]>(z)[0] is the real part of z and reinterpret_cast<T(&)[2]>(z)[1] is the imaginary part of z.
 * For any pointer to an element of an array of std::complex<T> named p and any valid array index i, reinterpret_cast<T*>(p)[2 * i] is the real part of the complex number p[i], and reinterpret_cast<T*>(p)[2 * i + 1] is the imaginary part of the complex number p[i].
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose imaginary component is extracted
 * @return imaginary component of argument
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline T& imag(std::complex<T>& z) {
  return reinterpret_cast<T(&)[2]>(z)[1];
}

/**
 * Return the imaginary component of the complex argument.
 * For any object z of type std::complex<T>, reinterpret_cast<T(&)[2]>(z)[0] is the real part of z and reinterpret_cast<T(&)[2]>(z)[1] is the imaginary part of z.
 * For any pointer to an element of an array of std::complex<T> named p and any valid array index i, reinterpret_cast<T*>(p)[2 * i] is the real part of the complex number p[i], and reinterpret_cast<T*>(p)[2 * i + 1] is the imaginary part of the complex number p[i].
 *
 * @tparam T value type of complex argument
 * @param[in] z complex value whose imaginary component is extracted
 * @return imaginary component of argument
 */
template <typename T, require_floating_point_t<T>* = nullptr>
inline T& imag(std::complex<T>&& z) {
  return reinterpret_cast<T(&)[2]>(z)[1];
}

}  // namespace math
}  // namespace stan

#endif
