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
 * @tparam Container type of the vector (must be derived from \c
 * Eigen::MatrixBase)
 * @param x Vector.
 * @return L1 norm of v.
 */
template <typename Container,
          require_eigen_vt<std::is_arithmetic, Container>* = nullptr>
inline double norm1(Container&& x) {
  return x.template lpNorm<1>();
}

/**
 * Returns L1 norm of a vector. For vectors that equals the
 * sum of magnitudes of its individual elements.
 *
 * @tparam Container type of the vector (must be derived from \c std::Vector)
 * @param x Vector.
 * @return L1 norm of v.
 */
template <typename Container, require_std_vector_t<Container>* = nullptr>
inline auto norm1(Container&& x) {
  return apply_vector_unary<Container>::reduce(
      std::forward<Container>(x),
      [](auto&& x_) { return norm1(std::forward<decltype(x_)>(x_)); });
}

}  // namespace math
}  // namespace stan

#endif
