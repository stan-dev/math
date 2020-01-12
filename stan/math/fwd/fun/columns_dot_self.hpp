#ifndef STAN_MATH_FWD_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_FWD_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/dot_self.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, 1, C> columns_dot_self(
    const Eigen::Matrix<fvar<T>, R, C>& x) {
  Eigen::Matrix<fvar<T>, 1, C> ret(1, x.cols());
  for (size_type i = 0; i < x.cols(); i++) {
    Eigen::Matrix<fvar<T>, R, 1> ccol = x.col(i);
    ret(0, i) = dot_self(ccol);
  }
  return ret;
}
}  // namespace math
}  // namespace stan
#endif
