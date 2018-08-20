#ifndef STAN_MATH_PRIM_SCAL_FUN_PROMOTE_DOUBLE_TO_T_HPP
#define STAN_MATH_PRIM_SCAL_FUN_PROMOTE_DOUBLE_TO_T_HPP

#include <vector>

namespace stan {
namespace math {

template <typename T, typename R>
inline const R& promote_double_to_T(const R& x) {
  return x;
}

template <typename T>
inline T promote_double_to_T(double x) {
  return x;
}

}  // namespace math
}  // namespace stan
#endif
