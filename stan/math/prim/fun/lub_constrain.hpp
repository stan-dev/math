#ifndef STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
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
template <typename T, typename L, typename U>
inline auto lub_constrain(T&& x, L&& lb, U&& ub) {
  const auto& x_ref = to_ref(x);
  const auto& lb_ref = to_ref(lb);
  const auto& ub_ref = to_ref(ub);
  check_less("lub_constrain", "lb", value_of(lb_ref), value_of(ub_ref));
  check_finite("lub_constrain", "lb", value_of(lb_ref));
  check_finite("lub_constrain", "ub", value_of(ub_ref));
  return eval(
      add(elt_multiply(subtract(ub_ref, lb_ref), inv_logit(x_ref)), lb_ref));
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
template <typename T, typename L, typename U>
inline auto lub_constrain(T&& x, L&& lb, U&& ub, return_type_t<T, L, U>& lp) {
  auto&& x_ref = to_ref(std::forward<T>(x));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  auto&& ub_ref = to_ref(std::forward<U>(ub));

  check_less("lub_constrain", "lb", value_of(lb_ref), value_of(ub_ref));
  check_finite("lub_constrain", "lb", value_of(lb_ref));
  check_finite("lub_constrain", "ub", value_of(ub_ref));

  auto diff = eval(subtract(std::forward<decltype(ub_ref)>(ub_ref), lb_ref));

  lp += sum(
      add(log(diff), subtract(-abs(x_ref), multiply(static_cast<double>(2),
                                                    log1p_exp(-abs(x_ref))))));

  return eval(add(elt_multiply(diff, inv_logit(x_ref)), lb_ref));
}

}  // namespace math
}  // namespace stan
#endif
