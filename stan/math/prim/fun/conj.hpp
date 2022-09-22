#ifndef STAN_MATH_PRIM_FUN_CONJ_HPP
#define STAN_MATH_PRIM_FUN_CONJ_HPP

#include <complex>
#include <stan/math/prim/meta.hpp>

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

/**
 * Return the complex conjugate the Eigen object.
 *
 * @tparam Eig A type derived from `Eigen::EigenBase`
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename Eig, require_eigen_vt<is_complex, Eig>* = nullptr>
inline auto conj(const Eig& z) {
  return z.conjugate();
}

/**
 * Return the complex conjugate the vector with complex scalar components.
 *
 * @tparam StdVec A `std::vector` type with complex scalar type
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename StdVec, require_std_vector_st<is_complex, StdVec>* = nullptr>
inline auto conj(const StdVec& z) {
  promote_scalar_t<scalar_type_t<StdVec>, StdVec> result(z.size());
  std::transform(z.begin(), z.end(), result.begin(),
                 [](auto&& x) { return stan::math::conj(x); });
  return result;
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
