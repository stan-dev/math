#ifndef STAN_MATH_PRIM_MAT_FUN_OFFSET_MULTIPLIER_CONSTRAIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_OFFSET_MULTIPLIER_CONSTRAIN_HPP

#include <stan/math/prim/mat/err/check_cholesky_factor.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_finite.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/diagonal.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * <p>If the offset is zero and the multiplier is one this
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of vector
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] x Unconstrained vector input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return linear transformed value correspdonding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S>
inline Eigen::Matrix<typename return_type<T, M, S>::type, -1, 1>
offset_multiplier_constrain(const Eigen::Matrix<T, -1, 1>& x,
                            const Eigen::Matrix<M, -1, 1>& mu,
                            const Eigen::Matrix<S, -1, -1>& sigma) {
  static const char* function = "offset_multiplier_constrain";
  check_finite(function, "offset", mu);
  check_finite(function, "multiplier", sigma);
  check_cholesky_factor(function, "multiplier", sigma);
  check_square(function, "multiplier", sigma);
  check_consistent_sizes(function, "multiplier", sigma.col(0),
                         "contrained vector", x);
  const size_t N = sigma.col(0).size();
  if (sigma == Eigen::Matrix<S, -1, -1>::Identity(N, N)) {
    if (mu == Eigen::Matrix<M, -1, 1>::Zero(N, 1))
      return x;
    return mu + x;
  }
  return mu + sigma * x;
}

/**
 * Return the linearly transformed value for the specified unconstrained input
 * and specified offset and multiplier, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = mu + sigma * x\f$
 *
 * <p>where mu is the offset and sigma is the multiplier.
 *
 * If the offset is zero and multiplier is one, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T type of vector
 * @tparam M type of offset
 * @tparam S type of multiplier
 * @param[in] x Unconstrained vector input
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @param[in,out] lp Reference to log probability to increment.
 * @return linear transformed value corresponding to inputs
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 */
template <typename T, typename M, typename S>
inline Eigen::Matrix<typename return_type<T, M, S>::type, -1, 1>
offset_multiplier_constrain(const Eigen::Matrix<T, -1, 1>& x,
                            const Eigen::Matrix<M, -1, 1>& mu,
                            const Eigen::Matrix<S, -1, -1>& sigma, T& lp) {
  static const char* function = "offset_multiplier_constrain";
  check_finite(function, "offset", mu);
  check_finite(function, "multiplier", sigma);
  check_cholesky_factor(function, "multiplier", sigma);
  check_square(function, "multiplier", sigma);
  check_consistent_sizes(function, "multiplier", sigma.col(0),
                         "contrained vector", x);
  const size_t N = sigma.col(0).size();
  if (sigma == Eigen::Matrix<S, -1, -1>::Identity(N, N)) {
    if (mu == Eigen::Matrix<M, -1, 1>::Zero(N, 1))
      return x;
    return mu + x;
  }
  lp += sum(log(diagonal(sigma)));
  return mu + sigma * x;
}

}  // namespace math

}  // namespace stan

#endif
