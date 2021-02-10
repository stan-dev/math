#ifndef STAN_MATH_PRIM_FUN_LB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unconstrained value that produces the specified
 * lower-bound constrained value.
 *
 * @tparam T type of bounded object
 * @tparam L type of lower bound
 * @param[in] y input object
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if y is lower than the lower bound
 */
template <typename T, typename L>
inline auto lb_free(T&& y, L&& lb) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  check_finite("lb_constrain", "lb", value_of(lb_ref));
  check_greater_or_equal("lb_free", "Lower bounded variable", value_of(y_ref),
                         value_of(lb_ref));
  return log(subtract(std::forward<decltype(y_ref)>(y_ref),
                      std::forward<decltype(lb_ref)>(lb_ref)));
}

}  // namespace math
}  // namespace stan
#endif
