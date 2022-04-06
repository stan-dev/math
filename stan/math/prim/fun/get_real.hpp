#ifndef STAN_MATH_PRIM_FUN_GET_REAL_HPP
#define STAN_MATH_PRIM_FUN_GET_REAL_HPP

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
template <typename T>
inline T get_real(const std::complex<T>& z) {
  return z.real();
}

/**
 * Return the real component of the complex argument.
 *
 * @tparam Eig A type derived from `Eigen::EigenBase`
 * @param[in] z complex value whose real component is extracted
 * @return real component of argument
 */
template <typename Eig, require_eigen_vt<is_complex, Eig>* = nullptr>
inline auto get_real(const Eig& z) {
  return z.real();
}

/**
 * Return the real component of the complex argument.
 *
 * @tparam StdVec A `std::vector` type with complex scalar type
 * @param[in] z complex value whose real component is extracted
 * @return real component of argument
 */
template <typename StdVec, require_std_vector_st<is_complex, StdVec>* = nullptr>
inline auto get_real(const StdVec& z) {
  promote_scalar_t<base_type_t<StdVec>, StdVec> result(z.size());
  std::transform(z.begin(), z.end(), result.begin(),
                 [](auto&& x) { return get_real(x); });
  return result;
}

}  // namespace math
}  // namespace stan

#endif
