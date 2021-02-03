#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb Lower bound
 * @return Constrained matrix
 */
template <typename T, typename L>
inline auto lb_constrain(T&& x, L&& lb) {
  return make_holder(
      [](const auto& lb_ref, const auto& xx) {
        check_finite("lb_constrain", "lb", value_of(lb_ref));
        return add(exp(xx), lb_ref);
      },
      to_ref(std::forward<L>(lb)), to_ref(std::forward<T>(x)));
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * @tparam T Type of Matrix
 * @tparam L type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained Matrix input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L>
inline auto lb_constrain(T&& x, L&& lb, return_type_t<T, L>& lp) {
  auto&& x_ref = to_ref(std::forward<T>(x));
  lp += sum(x_ref);
  return make_holder([](const auto& xx_ref, const auto& lb_ref) {
    check_finite("lb_constrain", "lb", value_of(lb_ref));
    return add(exp(xx_ref), lb_ref);
  }, std::forward<decltype(x_ref)>(x_ref),
  to_ref(std::forward<L>(lb)));
}

}  // namespace math
}  // namespace stan

#endif
