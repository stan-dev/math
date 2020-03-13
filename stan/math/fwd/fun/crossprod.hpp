#ifndef STAN_MATH_FWD_FUN_CROSSPROD_HPP
#define STAN_MATH_FWD_FUN_CROSSPROD_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/fwd/fun/multiply.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, C, C> crossprod(
    const Eigen::Matrix<fvar<T>, R, C>& m) {
  if (m.rows() == 0) {
    return {};
  }
  return multiply(transpose(m), m);
}

}  // namespace math
}  // namespace stan
#endif
