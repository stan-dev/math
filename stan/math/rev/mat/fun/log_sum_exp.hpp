#ifndef STAN_MATH_REV_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

namespace internal {

// these function and the following class just translate
// log_sum_exp for std::vector for Eigen::Matrix

template <int R, int C>
inline double log_sum_exp_as_double(const Eigen::Matrix<var, R, C>& x) {
  const double max = x.val().maxCoeff();
  if (!std::isfinite(max)) {
    return max;
  }
  return max + std::log((x.val().array() - max).exp().sum());
}

class log_sum_exp_matrix_vari : public op_matrix_vari {
 public:
  template <int R, int C>
  explicit log_sum_exp_matrix_vari(const Eigen::Matrix<var, R, C>& x)
      : op_matrix_vari(log_sum_exp_as_double(x), x) {}
  void chain() {
    Eigen::Map<vector_vi> vis_map(vis_, size_);
    vis_map.adj().array() += adj_ * (vis_map.val().array() - val_).exp();
  }
};
}  // namespace internal

/**
 * Returns the log sum of exponentials.
 *
 * @param x matrix
 */
template <int R, int C>
inline var log_sum_exp(const Eigen::Matrix<var, R, C>& x) {
  return var(new internal::log_sum_exp_matrix_vari(x));
}

}  // namespace math
}  // namespace stan
#endif
