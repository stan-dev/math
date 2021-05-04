#ifndef STAN_MATH_PRIM_SCAL_FUN_SIGNBIT_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SIGNBIT_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Return `true` if the specified argument is negative and `false`
 * otherwise.
 *
 * Overloads `std::signbit` from `<cmath>` for argument-dependent
 * lookup.
 *
 * @tparam ADType type of argument
 * @param[in] v argument
 * @return `true` if the argument is negative
 */
template <typename ADType, require_autodiff_t<ADType>* = nullptr>
inline bool signbit(ADType&& v) {
  using std::signbit;
  return signbit(v.val());
}

}  // namespace math
}  // namespace stan

#endif
