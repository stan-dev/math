#ifndef STAN_MATH_PRIM_SCAL_FUN_ISINF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISINF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return true if specified argument is infinite (positive or
 * negative).
 *
 * Overloads `std::isinf` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is infinite
 */
template <typename T, typename = require_autodiff_t<T>>
inline bool isinf(const T& v) {
  return is_inf(v);
}

}  // namespace math
}  // namespace stan

#endif
