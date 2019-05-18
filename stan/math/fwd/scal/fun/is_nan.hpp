#ifndef STAN_MATH_FWD_SCAL_FUN_IS_NAN_HPP
#define STAN_MATH_FWD_SCAL_FUN_IS_NAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Returns true if the input's value is NaN and false otherwise.
 *
 * Delegates to <code>is_nan</code>.
 *
 * @param x Value to test.
 * @return <code>true</code> if the value is NaN and <code>false</code> otherwise.
 */
template <typename T>
inline bool is_nan(const fvar<T>& x) {
  return is_nan(x.val());
}

}  // namespace math
}  // namespace stan
#endif
