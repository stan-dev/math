#ifndef STAN_MATH_PRIM_FUN_LUB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LUB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/lb_free.hpp>
#include <stan/math/prim/fun/ub_free.hpp>
#include <stan/math/prim/fun/subtract.hpp>

namespace stan {
namespace math {

/**
 * Return the unconstrained scalar that transforms to the
 * specified lower- and upper-bounded scalar given the specified
 * bounds.
 *
 * <p>The transform in <code>lub_constrain(T, double, double)</code>,
 * is reversed by a transformed and scaled logit,
 *
 * <p>\f$f^{-1}(y) = \mbox{logit}(\frac{y - L}{U - L})\f$
 *
 * where \f$U\f$ and \f$L\f$ are the lower and upper bounds.
 *
 * @tparam T type of bounded object
 * @tparam L type of lower bound
 * @tparam U type of upper bound
 * @param y constrained value
 * @param lb lower bound
 * @param ub upper bound
 * @return the free object that transforms to the input scalar
 *   given the bounds
 * @throw std::invalid_argument if the lower bound is greater than
 *   the upper bound, y is less than the lower bound, or y is
 *   greater than the upper bound
 */
template <typename T, typename L, typename U>
inline auto lub_free(T&& y, L&& lb, U&& ub) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  auto&& ub_ref = to_ref(std::forward<U>(ub));
  check_finite("lub_free", "lb", value_of(lb_ref));
  check_finite("lub_free", "ub", value_of(ub_ref));
  check_bounded("lub_free", "Bounded variable", value_of(y_ref),
                value_of(lb_ref), value_of(ub_ref));
  return eval(
      logit(divide(subtract(std::forward<decltype(y_ref)>(y_ref), lb_ref),
                   subtract(std::forward<decltype(ub_ref)>(ub_ref), lb_ref))));
}

}  // namespace math
}  // namespace stan
#endif
