#ifndef STAN_MATH_PRIM_SCAL_FUN_ISNAN_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISNAN_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return true if specified argument is not-a-number.
 *
 * Overloads `std::isnan` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return true if argument is not-a-number
 */
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isnan(ADType&& v) {
  using std::isnan;
  return isnan(v.val());
}

}  // namespace math
}  // namespace stan

#endif
