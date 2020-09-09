#ifndef STAN_MATH_PRIM_SCAL_FUN_ISNORMAL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISNORMAL_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return true if specified argument is normal.  A number is normal if
 * it is finite, non-zero and not subnormal.
 *
 * Overloads `std::isnormal` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is normal
 */
template <typename ADType, require_autodiff_t<ADType>* = nullptr>
inline bool isnormal(ADType&& v) {
  using std::isnormal;
  return isnormal(v.val());
}

}  // namespace math
}  // namespace stan

#endif
