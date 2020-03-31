#ifndef STAN_MATH_PRIM_FUN_LOGB_HPP
#define STAN_MATH_PRIM_FUN_LOGB_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the value of the unbiased radix-independent exponent from
 * the floating-point argument.
 *
 * Implementation note: This implementation works for all autodiff
 * types because it is a step function, so can return a `double` which
 * will produce the correct zero derivative if promoted to an autodiff
 * variable.
 *
 * @tparam T floating-point type
 * @param[in] x floating-point argument
 * @return unbiased radix-independent exponent of the argument
 */
template <typename T, typename = require_autodiff_t<T>>
double logb(const T& x) {
  return std::logb(value_of_rec(x));
}

}  // namespace math
}  // namespace stan

#endif
