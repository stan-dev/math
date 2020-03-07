#ifndef STAN_MATH_PRIM_SCAL_FUN_ISINF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ISINF_HPP

#include <stan/math/prim/meta.hpp>

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
template <typename ADType, require_autodiff_t<ADType>...>
inline bool isinf(ADType&& v) {
  using std::isinf;
  return isinf(v.val());
}

}  // namespace math
}  // namespace stan

#endif
