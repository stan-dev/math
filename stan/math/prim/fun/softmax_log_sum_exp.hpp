#ifndef STAN_MATH_PRIM_FUN_SOFTMAX_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_FUN_SOFTMAX_LOG_SUM_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/softmax.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Return the softmax and the log sum of exponentials of the specified
 * vector or matrix, using a single evaluation of the exponential per
 * element. For a matrix, both are taken over all elements.
 *
 * With \f$m = \max(x)\f$ and \f$\theta = \exp(x - m)\f$,
 *
 * \f$\mbox{softmax}(x) = \theta / \sum_{n=1}^N \theta_n\f$ and
 * \f$\mbox{log\_sum\_exp}(x) = m + \log \sum_{n=1}^N \theta_n\f$.
 *
 * The results match `softmax(x)` and `log_sum_exp(x)`, including for
 * non-finite input: if \f$m\f$ is not finite, the log sum of exponentials
 * is \f$m\f$ and the softmax is NaN.
 *
 * @tparam T type of Eigen vector or matrix
 * @param x Eigen vector or matrix
 * @return pair of the softmax of `x` and the log sum of exponentials of
 *   `x`. For empty input, returns an empty softmax and negative infinity
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline std::pair<plain_type_t<T>, double> softmax_log_sum_exp(const T& x) {
  if (x.size() == 0) {
    return {x, NEGATIVE_INFTY};
  }
  const auto& x_ref = to_ref(x);
  const double max = x_ref.maxCoeff();
  const auto theta = (x_ref.array() - max).exp().eval();
  const double theta_sum = theta.sum();
  const double lse = std::isfinite(max) ? max + std::log(theta_sum) : max;
  return {plain_type_t<T>((theta / theta_sum).matrix()), lse};
}

/**
 * Return the softmax and the log sum of exponentials of the specified
 * vector or matrix of autodiff values. This happens with nested autodiff,
 * e.g. the values of an `fvar<var>` vector, and falls back to the
 * autodiff specializations of `softmax()` and `log_sum_exp()`.
 *
 * @tparam T type of Eigen vector or matrix with autodiff scalar
 * @param x vector or matrix
 * @return pair of the softmax of `x` and the log sum of exponentials of `x`
 */
template <typename T, require_eigen_vt<is_autodiff, T>* = nullptr>
inline std::pair<plain_type_t<T>, value_type_t<T>> softmax_log_sum_exp(
    const T& x) {
  return {plain_type_t<T>(softmax(x)), log_sum_exp(x)};
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
