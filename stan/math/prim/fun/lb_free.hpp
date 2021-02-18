#ifndef STAN_MATH_PRIM_FUN_LB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
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
inline return_type_t<T, L> lb_free(const T& y, const L& lb) {
  using std::log;
  if (lb == NEGATIVE_INFTY) {
    return identity_free(y);
  }
  check_greater_or_equal("lb_free", "Lower bounded variable", y, lb);
  return log(y - lb);
}

}  // namespace math
}  // namespace stan
#endif
