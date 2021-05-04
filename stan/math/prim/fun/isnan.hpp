#ifndef STAN_MATH_PRIM_SCAL_FUN_ISNAN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISNAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] x argument
 * @return true if argument is not-a-number
 */
template <typename T, typename = require_autodiff_t<T>>
inline bool isnan(const T& x) {
  return is_nan(x);
}

}  // namespace math
}  // namespace stan

#endif
