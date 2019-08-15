#ifndef STAN_MATH_PRIM_MAT_FUN_ORDERED_FREE_HPP
#define STAN_MATH_PRIM_MAT_FUN_ORDERED_FREE_HPP


#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
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
 * @param y Vector of positive, ordered scalars.
 * @return Free vector that transforms into the input vector.
 * @tparam T Type of scalar.
 * @throw std::domain_error if y is not a vector of positive,
 *   ordered scalars.
 */
template <typename T, typename = enable_if_eigen<T>, std::enable_if_t<T::ColsAtCompileTime == 1>* = nullptr>
auto ordered_free(
    const T& y) {
  check_ordered("stan::math::ordered_free", "Ordered variable", y);
  auto k = y.size();
  typename T::PlainObject x(k);
  if (k == 0)
    return x;
  x[0] = y[0];
  for (int i = 1; i < k; ++i)
    x[i] = log(y[i] - y[i - 1]);
  return x;
}
}  // namespace math
}  // namespace stan
#endif
