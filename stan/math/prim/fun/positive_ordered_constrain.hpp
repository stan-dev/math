#ifndef STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_POSITIVE_ORDERED_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing positive ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> positive_ordered_constrain(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  using size_type = typename index_type<Matrix<T, Dynamic, 1>>::type;

  size_type k = x.size();
  Matrix<T, Dynamic, 1> y(k);
  if (k == 0) {
    return y;
  }
  y[0] = exp(x[0]);
  for (size_type i = 1; i < k; ++i) {
    y[i] = y[i - 1] + exp(x[i]);
  }
  return y;
}

/**
 * Return a positive valued, increasing positive ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @tparam T type of elements in the vector
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> positive_ordered_constrain(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, T& lp) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using size_type = typename index_type<Matrix<T, Dynamic, 1>>::type;

  for (size_type i = 0; i < x.size(); ++i) {
    lp += x(i);
  }
  return positive_ordered_constrain(x);
}

}  // namespace math
}  // namespace stan

#endif
