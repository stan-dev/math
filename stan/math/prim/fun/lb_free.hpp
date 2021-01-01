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
 * If the lower bound is negative infinity, it is ignored and
 * the function reduces to <code>identity_free(y)</code>.
 *
 * @tparam T type of scalar
 * @tparam L type of lower bound
 * @param[in] y input scalar
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if y is lower than the lower bound
 */
template <typename T, typename L>
inline auto lb_free(T&& y, const L& lb) {
  if (unlikely(is_negative_infinity(lb))) {
    return identity_free(std::forward<T>(y));
  }
  check_greater_or_equal("lb_free", "Lower bounded variable", value_of(y), value_of(lb));
  return log(subtract(y, lb));
}

}  // namespace math
}  // namespace stan
#endif
