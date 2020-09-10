#ifndef STAN_MATH_PRIM_SCAL_FUN_ISFINITE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISFINITE_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return true if specified argument is finite (not infinite and not
 * not-a-number).
 *
 * Overloads `std::isfinite` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is finite
 */
template <typename ADType, require_autodiff_t<ADType>* = nullptr>
inline bool isfinite(ADType&& v) {
  using std::isfinite;
  return isfinite(v.val());
}

}  // namespace math
}  // namespace stan

#endif
