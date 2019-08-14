#ifndef STAN_MATH_PRIM_MAT_FUN_TRANSPOSE_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRANSPOSE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

template <typename T, typename = enable_if_eigen<T>>
auto inline transpose(const T& m) {
  return m.transpose();
}

}  // namespace math
}  // namespace stan
#endif
