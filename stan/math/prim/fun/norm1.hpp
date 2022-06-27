#ifndef STAN_MATH_PRIM_FUN_NORM1_HPP
#define STAN_MATH_PRIM_FUN_NORM1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns L1 norm of a vector. For vectors that equals the
 * sum of magnitudes of its individual elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 * @return L1 norm of v.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline double norm1(const T& v) {
  ref_type_t<T> v_ref = v;
  return v_ref.template lpNorm<1>();
}

/**
 * Returns L1 norm of a vector. For vectors that equals the
 * sum of magnitudes of its individual elements.
 *
 * @tparam Container type of the vector (must be derived from \c std::Vector)
 * @param v Vector.
 * @return L1 norm of v.
 */
template <typename Container, require_std_vector_t<Container>* = nullptr>
inline auto norm1(const Container& x) {
  return apply_vector_unary<Container>::reduce(
      x, [](const auto& v) { return norm1(v); });
}

}  // namespace math
}  // namespace stan

#endif
