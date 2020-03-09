#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param[in] x free scalar.
 * @param[in] ub upper bound
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, require_stan_scalar_t<T>* = nullptr>
inline auto ub_constrain(T&& x, U&& ub) {
  using std::exp;
  if (ub == INFTY) {
    return identity_constrain(x, ub);
  }
  return ub - exp(x);
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, ub)</code>.
 *
 * @tparam EigT type derived from EigenBase
 * @tparam U type of upper bound
 * @param[in] x free scalar.
 * @param[in] ub upper bound
 * @return scalar constrained to have upper bound
 */
template <typename EigT, typename U, require_eigen_t<EigT>* = nullptr>
inline auto ub_constrain(EigT&& x, U&& ub) {
  return x.unaryExpr([&ub](auto&& x_iter){
    return ub_constrain(x_iter, ub);
  }).eval();
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, ub)</code>.
 *
 * @tparam Vec type of standard vector
 * @tparam U type of upper bound
 * @param[in] x free scalar.
 * @param[in] ub upper bound
 * @return scalar constrained to have upper bound
 */
template <typename Vec, typename U, require_std_vector_t<Vec>* = nullptr>
inline auto ub_constrain(Vec&& x, U&& ub) {
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(), [&ub](auto&& x_iter){
    return ub_constrain(x_iter, ub);
  });
  return ret_x;
}


/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, ub)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename U, typename S, require_stan_scalar_t<T>* = nullptr>
inline auto ub_constrain(T&& x, U&& ub, S& lp) {
  using std::exp;
  if (ub == INFTY) {
    return identity_constrain(x, ub);
  }
  lp += x;
  return ub - exp(x);
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, ub)</code>.
 *
 * @tparam EigT Type derived from EigenBase
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename EigT, typename U, typename S, require_eigen_t<EigT>* = nullptr>
inline auto ub_constrain(EigT&& x, U& ub, S& lp) {
  return x.unaryExpr([&ub, &lp](auto&& x_iter){
    return ub_constrain(x_iter, ub, lp);
  }).eval();
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, ub)</code>.
 *
 * @tparam Vec type of standard vector
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename Vec, typename U, typename S, require_std_vector_t<Vec>* = nullptr>
inline auto ub_constrain(Vec&& x, U&& ub, S& lp) {
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(), [&ub, &lp](auto&& x_iter){
    return ub_constrain(x_iter, ub, lp);
  });
  return ret_x;
}

}  // namespace math

}  // namespace stan

#endif
