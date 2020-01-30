#ifndef STAN_MATH_PRIM_FUN_TRANSPOSE_HPP
#define STAN_MATH_PRIM_FUN_TRANSPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T, typename = require_eigen_t<T>>
inline auto transpose(const T& m) {
  return m.transpose().eval();
}

}  // namespace math
}  // namespace stan
#endif
