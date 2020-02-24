#ifndef STAN_MATH_PRIM_FUN_DIMS_HPP
#define STAN_MATH_PRIM_FUN_DIMS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T, typename = require_stan_scalar_t<T>>
inline void dims(const T& x, std::vector<int>& result) {
  /* no op */
}
template <typename T, typename = require_eigen_t<T>, typename = void>
inline void dims(const T& x, std::vector<int>& result) {
  result.push_back(x.rows());
  result.push_back(x.cols());
}
template <typename T>
inline void dims(const std::vector<T>& x, std::vector<int>& result) {
  result.push_back(x.size());
  if (x.size() > 0) {
    dims(x[0], result);
  }
}

template <typename T>
inline std::vector<int> dims(const T& x) {
  std::vector<int> result;
  dims(x, result);
  return result;
}

}  // namespace math
}  // namespace stan
#endif
