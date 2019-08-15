#ifndef STAN_MATH_PRIM_MAT_FUN_ORDERED_CONSTRAIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_ORDERED_CONSTRAIN_HPP


#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an increasing ordered vector derived from the specified
 * free vector.  The returned constrained vector will have the
 * same dimensionality as the specified free vector.
 *
 * @param x Free vector of scalars.
 * @return Positive, increasing ordered vector.
 * @tparam T Type of scalar.
 */
template <typename T, typename = enable_if_eigen<T>, std::enable_if_t<T::ColsAtCompileTime == 1>* = nullptr>
auto ordered_constrain(
    const T& x) {
  using std::exp;

  auto k = x.size();
  auto y = x.array().exp().eval();
  y[0] = x[0];
  for (int i = 1; i < k; i++) {
    y[i] += y[i - 1];
  }
  return y;
}

/**
 * Return a positive valued, increasing ordered vector derived
 * from the specified free vector and increment the specified log
 * probability reference with the log absolute Jacobian determinant
 * of the transform.  The returned constrained vector
 * will have the same dimensionality as the specified free vector.
 *
 * @param x Free vector of scalars.
 * @param lp Log probability reference.
 * @return Positive, increasing ordered vector.
 * @tparam T Type of scalar.
 */
 template <typename T1, typename T2, typename = enable_if_eigen<T1>, std::enable_if_t<T1::ColsAtCompileTime == 1>* = nullptr>
inline auto ordered_constrain(
    const T1& x, T2& lp) {
  for (int i = 1; i < x.size(); ++i) {
    lp += x(i);
  }
  return ordered_constrain(x);
}

}  // namespace math

}  // namespace stan

#endif
