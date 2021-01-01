#ifndef STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/ub_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower- and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds.
 *
 * <p>The transform is the transformed and scaled inverse logit,
 *
 * <p>\f$f(x) = L + (U - L) \mbox{logit}^{-1}(x)\f$
 *
 * If the lower bound is negative infinity and upper bound finite,
 * this function reduces to <code>ub_constrain(x, ub)</code>.  If
 * the upper bound is positive infinity and the lower bound
 * finite, this function reduces to
 * <code>lb_constrain(x, lb)</code>.  If the upper bound is
 * positive infinity and the lower bound negative infinity,
 * this function reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T Type of scalar.
 * @tparam L Type of lower bound.
 * @tparam U Type of upper bound.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, L, U> lub_constrain(const T& x, const L& lb,
                                            const U& ub) {
  check_less("lub_constrain", "lb", lb, ub);
  if (lb == NEGATIVE_INFTY) {
    return ub_constrain(x, ub);
  }
  if (ub == INFTY) {
    return lb_constrain(x, lb);
  }
  return fma(ub - lb, inv_logit(x), lb);
}

/**
 * Return the lower- and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds.
 *
 * <p>The transform is the transformed and scaled inverse logit,
 *
 * <p>\f$f(x) = L + (U - L) \mbox{logit}^{-1}(x)\f$
 *
 * If the lower bound is negative infinity and upper bound finite,
 * this function reduces to <code>ub_constrain(x, ub)</code>.  If
 * the upper bound is positive infinity and the lower bound
 * finite, this function reduces to
 * <code>lb_constrain(x, lb)</code>.  If the upper bound is
 * positive infinity and the lower bound negative infinity,
 * this function reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T Type of scalar.
 * @tparam L Type of lower bound.
 * @tparam U Type of upper bound.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U, require_matrix_t<T>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb, const U& ub) {
  check_less("lub_constrain", "lb", lb, ub);
  if (unlikely(is_negative_infinity(lb))) {
    return ub_constrain(x, ub);
  }
  if (unlikely(is_infinity(ub))) {
    return lb_constrain(x, lb);
  }
  return add(subtract(ub, lb), elt_multiply(inv_logit(x), lb)).eval();
}

/**
 * Return the lower- and upper-bounded scalar derived by
 * transforming the specified free scalar given the specified
 * lower and upper bounds and increment the specified log
 * density with the log absolute Jacobian determinant.
 *
 * <p>The transform is as defined in
 * <code>lub_constrain(T, double, double)</code>.  The log absolute
 * Jacobian determinant is given by
 *
 * <p>\f$\log \left| \frac{d}{dx} \left(
 *                L + (U-L) \mbox{logit}^{-1}(x) \right)
 *            \right|\f$
 *
 * <p>\f$ {} = \log |
 *         (U-L)
 *         \, (\mbox{logit}^{-1}(x))
 *         \, (1 - \mbox{logit}^{-1}(x)) |\f$
 *
 * <p>\f$ {} = \log (U - L) + \log (\mbox{logit}^{-1}(x))
 *                          + \log (1 - \mbox{logit}^{-1}(x))\f$
 *
 * <p>If the lower bound is negative infinity and upper bound finite,
 * this function reduces to <code>ub_constrain(x, ub, lp)</code>.  If
 * the upper bound is positive infinity and the lower bound
 * finite, this function reduces to
 * <code>lb_constrain(x, lb, lp)</code>.  If the upper bound is
 * positive infinity and the lower bound negative infinity,
 * this function reduces to <code>identity_constrain(x, lp)</code>.
 *
 * @tparam T Type of scalar.
 * @tparam L Type of lower bound.
 * @tparam U Type of upper bound.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @param[in,out] lp Log probability scalar reference.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U, typename S, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T, L, U> lub_constrain(const T& x, const L& lb,
                                            const U& ub, S& lp) {
  using std::exp;
  using std::log;
  check_less("lub_constrain", "lb", lb, ub);
  if (lb == NEGATIVE_INFTY) {
    return ub_constrain(x, ub, lp);
  }
  if (ub == INFTY) {
    return lb_constrain(x, lb, lp);
  }
  if (x > 0) {
    lp += log(ub - lb) - x - 2 * log1p(exp(-x));
  } else {
    lp += log(ub - lb) + x - 2 * log1p(exp(x));
  }
  return fma(ub - lb, inv_logit(x), lb);
}

template <typename T, typename L, typename U, typename S, require_matrix_t<T>* = nullptr>
inline auto lub_constrain(const T& x, const L& lb,
                                            const U& ub, S& lp) {
  check_less("lub_constrain", "lb", lb, ub);
  if (unlikely(is_negative_infinity(lb))) {
    return ub_constrain(x, ub, lp);
  }
  if (unlikely(is_infinity(ub))) {
    return lb_constrain(x, lb, lp);
  }
  auto x_neg_abs = -abs(x);
  lp += sum(subtract(add(log(subtract(ub, lb)), x), multiply(2, log1p(exp(x)))));
  return add(subtract(ub, lb), elt_multiply(inv_logit(x), lb)).eval();
}

}  // namespace math
}  // namespace stan
#endif
