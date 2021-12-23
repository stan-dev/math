#ifndef STAN_MATH_PRIM_FUN_CONJ_HPP
#define STAN_MATH_PRIM_FUN_CONJ_HPP

#include <complex>

namespace stan {
namespace math {

/**
 * Return the complex conjugate the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename V>
inline std::complex<V> conj(const std::complex<V>& z) {
  return std::conj(z);
}

namespace internal {
/**
 * Return the complex conjugate the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename V>
inline std::complex<V> complex_conj(const std::complex<V>& z) {
  return {z.real(), -z.imag()};
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
