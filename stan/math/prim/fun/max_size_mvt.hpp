#ifndef STAN_MATH_PRIM_FUN_MAX_SIZE_MVT_HPP
#define STAN_MATH_PRIM_FUN_MAX_SIZE_MVT_HPP

#include <stan/math/prim/fun/size_mvt.hpp>
#include <cstdlib>

namespace stan {
namespace math {

template <typename T1, typename T2>
size_t max_size_mvt(const T1& x1, const T2& x2) {
  size_t result = size_mvt(x1);
  result = result > size_mvt(x2) ? result : size_mvt(x2);
  return result;
}

template <typename T1, typename T2, typename T3>
size_t max_size_mvt(const T1& x1, const T2& x2, const T3& x3) {
  size_t result = size_mvt(x1);
  result = result > size_mvt(x2) ? result : size_mvt(x2);
  result = result > size_mvt(x3) ? result : size_mvt(x3);
  return result;
}

template <typename T1, typename T2, typename T3, typename T4>
size_t max_size_mvt(const T1& x1, const T2& x2, const T3& x3, const T4& x4) {
  size_t result = size_mvt(x1);
  result = result > size_mvt(x2) ? result : size_mvt(x2);
  result = result > size_mvt(x3) ? result : size_mvt(x3);
  result = result > size_mvt(x4) ? result : size_mvt(x4);
  return result;
}

}  // namespace math
}  // namespace stan
#endif
