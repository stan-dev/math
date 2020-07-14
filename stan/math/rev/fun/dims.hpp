#ifndef STAN_MATH_REV_FUN_DIMS_HPP
#define STAN_MATH_REV_FUN_DIMS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T, typename = require_eigen_t<T>, typename = void>
inline void dims(const var_value<T>& x, std::vector<int>& result) {
  dims(*x.vi_, result);
}
template <typename T, typename = require_eigen_t<T>, typename = void>
inline void dims(const vari_value<T>& x, std::vector<int>& result) {
  result.push_back(x.rows());
  result.push_back(x.cols());
}

}  // namespace math
}  // namespace stan

#endif  // DIMS_HPP
