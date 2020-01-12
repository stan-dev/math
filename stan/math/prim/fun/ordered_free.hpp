#ifndef STAN_MATH_PRIM_FUN_ORDERED_FREE_HPP
#define STAN_MATH_PRIM_FUN_ORDERED_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the vector of unconstrained scalars that transform to
 * the specified positive ordered vector.
 *
 * <p>This function inverts the constraining operation defined in
 * <code>ordered_constrain(Matrix)</code>,
 *
 * @tparam T type of elements in the vector
 * @param y Vector of positive, ordered scalars.
 * @return Free vector that transforms into the input vector.
 * @throw std::domain_error if y is not a vector of positive,
 *   ordered scalars.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> ordered_free(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& y) {
  check_ordered("stan::math::ordered_free", "Ordered variable", y);
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::log;
  using size_type = typename index_type<Matrix<T, Dynamic, 1>>::type;

  size_type k = y.size();
  Matrix<T, Dynamic, 1> x(k);
  if (k == 0) {
    return x;
  }
  x[0] = y[0];
  for (size_type i = 1; i < k; ++i) {
    x[i] = log(y[i] - y[i - 1]);
  }
  return x;
}

}  // namespace math
}  // namespace stan

#endif
