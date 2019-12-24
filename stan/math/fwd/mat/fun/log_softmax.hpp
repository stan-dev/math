#ifndef STAN_MATH_FWD_MAT_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_FWD_MAT_FUN_LOG_SOFTMAX_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/mat/fun/softmax.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/log_softmax.hpp>
#include <stan/math/prim/mat/fun/softmax.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

template <typename T, require_t<is_fvar<scalar_type_t<T>>>...>
inline auto log_softmax(T&& x) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [&](auto& alpha){
    using Eigen::Dynamic;
    using Eigen::Matrix;

    using T_value_type = value_type_t<T>;
    using T_scalar_type = typename T_value_type::Scalar;

    Matrix<T_scalar_type, Dynamic, 1> alpha_t = alpha.val();
    Matrix<T_scalar_type, Dynamic, 1> softmax_alpha_t = softmax(alpha_t);

    Matrix<T_value_type, Dynamic, 1> log_softmax_alpha(alpha.size());
    log_softmax_alpha.val() = log_softmax(alpha_t);
    log_softmax_alpha.d().setZero();

    for (int m = 0; m < alpha.size(); ++m) {
      T_scalar_type negative_alpha_m_d_times_softmax_alpha_t_m
          = -alpha(m).d_ * softmax_alpha_t(m);
      for (int k = 0; k < alpha.size(); ++k) {
        if (m == k) {
          log_softmax_alpha(k).d_
              += alpha(m).d_ + negative_alpha_m_d_times_softmax_alpha_t_m;
        } else {
          log_softmax_alpha(k).d_ += negative_alpha_m_d_times_softmax_alpha_t_m;
        }
      }
    }

    return log_softmax_alpha;
  });
}

}  // namespace math
}  // namespace stan
#endif
