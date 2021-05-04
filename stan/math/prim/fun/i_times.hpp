#ifndef STAN_MATH_PRIM_SCAL_FUN_I_TIMES_HPP
#define STAN_MATH_PRIM_SCAL_FUN_I_TIMES_HPP

#include <complex>

namespace stan {
namespace math {

/**
 * Return the specified complex number multiplied by `i`.
 *
 * This compound function is more efficient than mulitplying by a
 * constant `i` because it involves only a single arithmetic negation.
 *
 * @tparam value type of complex argument
 * @param[in] z complex argument
 * @return argument multipled by `i`
 */
template <typename T>
inline std::complex<T> i_times(const std::complex<T>& z) {
  return {-z.imag(), z.real()};
}

/**
 * Return the specified complex number multiplied by `-i`.
 *
 * This compound function is more efficient than mulitplying by the
 * constant `-i` because it involves only a single arithmetic
 * negation.
 *
 * @tparam value type of complex argument
 * @param[in] z complex argument
 * @return argument multipled by `-i`
 */
template <typename T>
inline std::complex<T> neg_i_times(const std::complex<T>& z) {
  return {z.imag(), -z.real()};
}

/**
 * Return the complex negation of the specified complex argument.
 *
 * @tparam V value type of complex argument
 * @param[in] z argument
 * @return negation of argument
 */
template <typename V>
inline std::complex<V> complex_negate(const std::complex<V>& z) {
  return {-z.real(), -z.imag()};
}

}  // namespace math
}  // namespace stan

#endif
