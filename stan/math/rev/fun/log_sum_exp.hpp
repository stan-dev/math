#ifndef STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
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
    avi_->adj_ += adj_ * inv_logit(avi_->val_ - bvi_->val_);
    bvi_->adj_ += adj_ * inv_logit(bvi_->val_ - avi_->val_);
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
      avi_->adj_ += adj_ * inv_logit(avi_->val_ - bd_);
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

/**
 * Returns the log sum of exponentials of the input.
 *
 * @tparam T A type inheriting from EigenBase with scalar type var
 * @param v input
 */
template <typename T, require_eigen_st<is_var, T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline var log_sum_exp(const T& v) {
  arena_t<decltype(v)> arena_v = v;
  arena_t<decltype(v.val())> arena_v_val = arena_v.val();
  var res = log_sum_exp(arena_v_val);

  reverse_pass_callback([arena_v, arena_v_val, res]() mutable {
    arena_v.adj()
        += res.adj() * (arena_v_val.array().val() - res.val()).exp().matrix();
  });

  return res;
}

/**
 * Returns the log sum of exponentials of the input.
 *
 * @tparam T A `var_value` with an input vector or matrix
 * @param x input
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var log_sum_exp(const T& x) {
  return make_callback_vari(log_sum_exp(x.val()), [x](const auto& res) mutable {
    x.adj() += res.adj() * (x.val().array().val() - res.val()).exp().matrix();
  });
}

/**
 * Returns the log sum of exponentials.
 *
 * @tparam T Type of input vector or matrix.
 * @param x matrix
 */
template <typename T, require_std_vector_st<is_var, T>* = nullptr>
inline auto log_sum_exp(const T& x) {
  return apply_vector_unary<T>::reduce(
      x, [](const auto& v) { return log_sum_exp(v); });
}

}  // namespace math
}  // namespace stan
#endif
