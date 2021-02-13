#ifndef STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LUB_CONSTRAIN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lb_constrain.hpp>
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
 * @tparam T A scalar or Eigen matrix.
 * @tparam L A scalar or Eigen matrix.
 * @tparam U A scalar or Eigen matrix.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U>* = nullptr, require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(T&& x, L&& lb, U&& ub) {
  check_less("lub_constrain", "lb", value_of(lb), value_of(ub));
  if (!is_positive_infinity(ub)) {
    return eval(lb_constrain(identity_constrain(x, ub), lb));
  } else if (!is_negative_infinity(lb)) {
    return eval(ub_constrain(identity_constrain(x, lb), ub));
  } else {
    return eval(
        add(elt_multiply(subtract(ub, lb), inv_logit(x)), lb));
  }
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
 * @tparam T A scalar or Eigen matrix.
 * @tparam L A scalar or Eigen matrix.
 * @tparam U A scalar or Eigen matrix.
 * @param[in] x Free scalar to transform.
 * @param[in] lb Lower bound.
 * @param[in] ub Upper bound.
 * @param[in,out] lp Log probability scalar reference.
 * @return Lower- and upper-bounded scalar derived from transforming
 *   the free scalar.
 * @throw std::domain_error if ub <= lb
 */
template <typename T, typename L, typename U, require_all_stan_scalar_t<T, L, U>* = nullptr, require_not_var_t<return_type_t<T, L, U>>* = nullptr>
inline auto lub_constrain(T&& x, L&& lb, U&& ub, return_type_t<T, L, U>& lp) {
  check_less("lub_constrain", "lb", value_of(lb), value_of(ub));
  if (!is_positive_infinity(ub)) {
    return lb_constrain(identity_constrain(x, ub), lb, lp);
  } else if (!is_negative_infinity(lb)) {
    return lb_constrain(identity_constrain(x, ub), lb, lp);
  } else {
    auto diff = subtract(std::forward<decltype(ub)>(ub), lb);
    lp += sum(
        add(log(diff), subtract(-abs(x), multiply(static_cast<double>(2),
                                                      log1p_exp(-abs(x))))));

    return add(elt_multiply(diff, inv_logit(x)), lb);
  }
}



}  // namespace math
}  // namespace stan
#endif
