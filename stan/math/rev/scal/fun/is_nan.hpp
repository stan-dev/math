#ifndef STAN_MATH_REV_SCAL_FUN_IS_NAN_HPP
#define STAN_MATH_REV_SCAL_FUN_IS_NAN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Returns true if the input's value is NaN and false otherwise.
 *
 * Delegates to <code>is_nan(double)</code>.
 *
 * @param x Value to test.
 *
 * @return <code>true</code> if the value is NaN and <code>false</code>
 * otherwise.
 */
template <typename T>
inline typename std::enable_if_t<is_var<T>::value, bool> is_nan(const T& x) {
  return is_nan(x.val());
}

}  // namespace math
}  // namespace stan
#endif
