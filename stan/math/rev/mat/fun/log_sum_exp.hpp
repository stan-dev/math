#ifndef STAN_MATH_REV_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_REV_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/vectorize/apply_vector_unary.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

namespace internal {

class log_sum_exp_matrix_vari : public op_matrix_vari {
 public:
  template <typename T>
  explicit log_sum_exp_matrix_vari(T&& x)
      : op_matrix_vari(log_sum_exp(std::forward<T>(x).val()),
                       std::forward<T>(x)) {}
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
inline auto log_sum_exp(T&& x) {
  return apply_vector_unary<T>::reduce(std::forward<T>(x), [&](auto& v) {
    return var(new internal::log_sum_exp_matrix_vari(v));
  });
}

}  // namespace math
}  // namespace stan
#endif
