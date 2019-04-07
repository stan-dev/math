#ifndef STAN_MATH_FWD_MAT_FUN_DETERMINANT_HPP
#define STAN_MATH_FWD_MAT_FUN_DETERMINANT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline fvar<T> determinant(const Eigen::Matrix<fvar<T>, R, C>& m) {
  check_square("determinant", "m", m);

  fvar<T> result;
  result.val_ = m.val_().determinant();
  result.d_ = result.val_ * (m.val_().inverse() * m.d_()).trace();

  return result;
}

}  // namespace math
}  // namespace stan
#endif
