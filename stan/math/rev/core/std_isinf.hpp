#ifndef STAN_MATH_REV_CORE_STD_ISINF_HPP
#define STAN_MATH_REV_CORE_STD_ISINF_HPP

#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/rev/core/std_numeric_limits.hpp>
#include <stan/math/rev/core/var.hpp>

namespace std {

/**
 * Return <code>1</code> if the specified argument is positive
 * infinity or negative infinity and <code>0</code> otherwise.
 *
 * @param a Argument.
 * @return 1 if argument is infinite and 0 otherwise.
 */
inline int isinf(const stan::math::var& a) {
  return stan::math::is_inf(a.val());
}

}  // namespace std

namespace stan {
namespace math {

// forwarding for ADL
inline auto isinf(const var& a) { return std::isinf(a); }

}  // namespace math
}  // namespace stan
#endif
