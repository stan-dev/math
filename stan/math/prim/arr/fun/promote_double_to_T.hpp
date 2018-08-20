#ifndef STAN_MATH_PRIM_ARR_FUN_PROMOTE_DOUBLE_TO_T_HPP
#define STAN_MATH_PRIM_ARR_FUN_PROMOTE_DOUBLE_TO_T_HPP

#include <vector>

namespace stan {
namespace math {

template <typename T, typename R>
inline const std::vector<R>& promote_double_to_T(const std::vector<R>& x) {
  return x;
}

template <typename T>
inline std::vector<T> promote_double_to_T(const std::vector<double>& x) {
  std::vector<T> out;
  out.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    out.push_back(x[i]);
  }
  return out;
}

}  // namespace math
}  // namespace stan
#endif
