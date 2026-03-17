#ifndef STAN_MATH_PRIM_FUN_NORM2_HPP
#define STAN_MATH_PRIM_FUN_NORM2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam Container type of the vector (must be derived from \c
 * Eigen::MatrixBase)
 * @param x Vector.
 * @return L2 norm of x.
 */
template <typename Container,
          require_eigen_vt<std::is_arithmetic, Container>* = nullptr>
inline double norm2(Container&& x) {
  return x.template lpNorm<2>();
}

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam Container type of the vector (must be derived from \c std::vector)
 * @param x Vector.
 * @return L2 norm of x.
 */
template <typename Container, require_std_vector_t<Container>* = nullptr>
inline auto norm2(Container&& x) {
  return apply_vector_unary<Container>::reduce(
      std::forward<Container>(x),
      [](auto&& x_) { return norm2(std::forward<decltype(x_)>(x_)); });
}

}  // namespace math
}  // namespace stan

#endif
