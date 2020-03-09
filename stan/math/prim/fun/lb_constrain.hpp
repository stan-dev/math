#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + U\f$
 *
 * <p>where \f$U\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of lower bound
 * @param[in] x Unconstrained scalar input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename U, require_stan_scalar_t<T>* = nullptr>
inline auto lb_constrain(T&& x, U&& lb) {
  using std::exp;
  if (lb == NEGATIVE_INFTY) {
    return identity_constrain(x, lb);
  }
  return exp(x) + lb;
}

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + U\f$
 *
 * <p>where \f$U\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam EigT type derived from EigenBase
 * @tparam U type of lower bound
 * @param[in] x Unconstrained scalar input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename EigT, typename U, require_eigen_t<EigT>* = nullptr>
inline auto lb_constrain(EigT&& x, U&& lb) {
  return x.unaryExpr([&lb](auto&& x_iter) { return lb_constrain(x_iter, lb); })
      .eval();
}

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + U\f$
 *
 * <p>where \f$U\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam Vec type of standard vector
 * @tparam U type of lower bound
 * @param[in] x Unconstrained scalar input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename Vec, typename U, require_std_vector_t<Vec>* = nullptr>
inline auto lb_constrain(Vec&& x, U&& lb) {
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(),
                 [&lb](auto&& x_iter) { return lb_constrain(x_iter, lb); });
  return ret_x;
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained scalar input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename U, typename S,
          require_stan_scalar_t<T>* = nullptr>
inline auto lb_constrain(T&& x, U&& lb, S& lp) {
  using std::exp;
  if (lb == NEGATIVE_INFTY) {
    return identity_constrain(x, lb);
  }
  lp += x;
  return exp(x) + lb;
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam EigT type derived from EigenBase
 * @tparam U type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained scalar input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename EigT, typename U, typename S,
          require_eigen_t<EigT>* = nullptr>
inline auto lb_constrain(EigT&& x, U&& lb, S& lp) {
  return x
      .unaryExpr(
          [&lb, &lp](auto&& x_iter) { return lb_constrain(x_iter, lb, lp); })
      .eval();
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam Vec standard vector type
 * @tparam U type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained scalar input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename Vec, typename U, typename S,
          require_std_vector_t<Vec>* = nullptr>
inline auto lb_constrain(Vec&& x, U&& lb, S& lp) {
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(), [&lb, &lp](auto&& x_iter) {
    return lb_constrain(x_iter, lb, lp);
  });
  return ret_x;
}

}  // namespace math

}  // namespace stan

#endif
