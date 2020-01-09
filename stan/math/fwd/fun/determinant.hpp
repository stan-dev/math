#ifndef STAN_MATH_FWD_FUN_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_DETERMINANT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline fvar<T> determinant(const Eigen::Matrix<fvar<T>, R, C>& m) {
  check_square("determinant", "m", m);

  const T vals = m.val().determinant();
  return fvar<T>(vals, vals * (m.val().inverse() * m.d()).trace());
}

}  // namespace math
}  // namespace stan
#endif
