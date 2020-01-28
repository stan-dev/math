#ifndef STAN_MATH_PRIM_FUN_REVERSE_HPP
#define STAN_MATH_PRIM_FUN_REVERSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/vectorize/apply_vector_unary.hpp>

namespace stan {
namespace math {

/**
 * Return a copy of the specified vector, row vector or array
 * in reversed order.
 *
 * This implicitly allows for reverse with nested arrays (e.g.
 * reverse(std::vector<VectorXd>)) which returns an array with each vector
 * reversed.
 *
 * @tparam T type of container
 * @param x container to reverse
 * @return Container in reversed order.
 */
template <typename T>
inline auto reverse(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.reverse(); });
}

}  // namespace math
}  // namespace stan
#endif
