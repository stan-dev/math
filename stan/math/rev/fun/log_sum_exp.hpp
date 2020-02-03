#ifndef STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/calculate_chain.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

namespace internal {

class log_sum_exp_vv_vari : public op_vv_vari {
 public:
  log_sum_exp_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(log_sum_exp(avi->val_, bvi->val_), avi, bvi) {}
  void chain() {
    avi_->adj_ += adj_ * calculate_chain(avi_->val_, val_);
    bvi_->adj_ += adj_ * calculate_chain(bvi_->val_, val_);
  }
};
class log_sum_exp_vd_vari : public op_vd_vari {
 public:
  log_sum_exp_vd_vari(vari* avi, double b)
      : op_vd_vari(log_sum_exp(avi->val_, b), avi, b) {}
  void chain() {
    if (val_ == NEGATIVE_INFTY) {
      avi_->adj_ += adj_;
    } else {
      avi_->adj_ += adj_ * calculate_chain(avi_->val_, val_);
    }
  }
};

}  // namespace internal

/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(const var& a, const var& b) {
  return var(new internal::log_sum_exp_vv_vari(a.vi_, b.vi_));
}
/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(const var& a, double b) {
  return var(new internal::log_sum_exp_vd_vari(a.vi_, b));
}
/**
 * Returns the log sum of exponentials.
 */
inline var log_sum_exp(double a, const var& b) {
  return var(new internal::log_sum_exp_vd_vari(b.vi_, a));
}

namespace internal {

class log_sum_exp_matrix_vari : public op_matrix_vari {
 public:
  template <typename T>
  explicit log_sum_exp_matrix_vari(const T& x)
      : op_matrix_vari(log_sum_exp(x.val()), x) {}
  void chain() {
    Eigen::Map<vector_vi> vis_map(vis_, size_);
    vis_map.adj().array() += adj_ * (vis_map.val().array() - val_).exp();
  }
};
}  // namespace internal

/**
 * Returns the log sum of exponentials.
 *
 * @tparam T Type of input vector or matrix.
 * @param x matrix
 */
template <typename T, require_t<is_var<scalar_type_t<T>>>...>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<T>::reduce(x, [&](const auto& v) {
    return var(new internal::log_sum_exp_matrix_vari(v));
  });
}

}  // namespace math
}  // namespace stan
#endif
