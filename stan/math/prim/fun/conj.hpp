#ifndef STAN_MATH_PRIM_FUN_CONJ_HPP
#define STAN_MATH_PRIM_FUN_CONJ_HPP

#include <stan/math/prim/meta.hpp>
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
template <typename V, require_complex_bt<std::is_arithmetic, V>* = nullptr>
inline auto conj(const V& z) {
  return std::conj(z);
}

/**
 * Return the complex conjugate the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename V, require_complex_bt<is_autodiff, V>* = nullptr>
inline V conj(const V& z) {
  return {z.real(), -z.imag()};
}

/**
 * Return the complex conjugate the Eigen object.
 *
 * @tparam Eig A type derived from `Eigen::EigenBase` with an inner complex
 * autodiff type.
 * @param[in] z argument
 * @return complex conjugate of the argument
 */
template <typename Eig, require_eigen_vt<is_complex, Eig>* = nullptr>
inline auto conj(Eig&& z) {
  return make_holder(
      [](auto&& z_inner) {
        return z_inner.unaryExpr([](auto&& z_i) { return conj(z_i); });
      },
      std::forward<Eig>(z));
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
  const auto z_size = z.size();
  promote_scalar_t<scalar_type_t<StdVec>, StdVec> result(z_size);
  for (std::size_t i = 0; i < z_size; ++i) {
    result[i] = conj(z[i]);
  }
  return result;
}

}  // namespace math
}  // namespace stan

#endif
