#ifndef STAN_MATH_REV_FUN_IS_NAN_HPP
#define STAN_MATH_REV_FUN_IS_NAN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Returns 1 if the input's value is NaN and 0 otherwise.
 *
 * Delegates to <code>is_nan(double)</code>.
 *
 * @tparam T type of input
 * @param v value to test
 * @return <code>1</code> if the value is NaN and <code>0</code> otherwise.
 */
inline bool is_nan(const var& v) { return is_nan(v.val()); }

}  // namespace math
}  // namespace stan
#endif
