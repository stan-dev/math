#ifndef STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_FWD_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
fvar<T> log_sum_exp(const Eigen::Matrix<fvar<T>, R, C>& v) {
  Eigen::Matrix<T, R, C> vals = v.val();
  Eigen::Matrix<T, R, C> exp_vals = vals.array().exp();

  return fvar<T>(log_sum_exp(vals),
                 v.d().cwiseProduct(exp_vals).sum() / exp_vals.sum());
}

}  // namespace math
}  // namespace stan
#endif
