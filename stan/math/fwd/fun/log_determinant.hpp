#ifndef STAN_MATH_FWD_FUN_LOG_DETERMINANT_HPP
#define STAN_MATH_FWD_FUN_LOG_DETERMINANT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/determinant.hpp>
#include <stan/math/fwd/fun/fabs.hpp>
#include <stan/math/fwd/fun/log.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline fvar<T> log_determinant(const Eigen::Matrix<fvar<T>, R, C>& m) {
  check_square("log_determinant", "m", m);

  return log(fabs(determinant(m)));
}

}  // namespace math
}  // namespace stan
#endif
