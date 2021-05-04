#ifndef STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_LOG_SOFTMAX_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/softmax.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log_softmax.hpp>
#include <stan/math/prim/fun/softmax.hpp>

namespace stan {
namespace math {

/**
 * Return the log softmax of the specified vector or container of vectors.
 *
 * @tparam T Type of input vector or matrix.
 * @param[in] x Unconstrained input vector.
 * @return Softmax of the input.
 * @throw std::domain_error If the input vector is size 0.
 */
template <typename T, require_vector_st<is_fvar, T>* = nullptr>
inline auto log_softmax(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& alpha) {
    using T_alpha = decltype(alpha);
    using T_fvar = value_type_t<T_alpha>;
    using T_fvar_inner = typename T_fvar::Scalar;

    const Eigen::Ref<const plain_type_t<T_alpha>>& alpha_ref = alpha;
    Eigen::Matrix<T_fvar_inner, -1, 1> alpha_t = alpha_ref.val();
    Eigen::Matrix<T_fvar_inner, -1, 1> softmax_alpha_t = softmax(alpha_t);

    Eigen::Matrix<T_fvar, -1, 1> log_softmax_alpha(alpha.size());
    log_softmax_alpha.val() = log_softmax(alpha_t);
    log_softmax_alpha.d().setZero();

    for (int m = 0; m < alpha.size(); ++m) {
      T_fvar_inner negative_alpha_m_d_times_softmax_alpha_t_m
          = -alpha_ref.coeff(m).d_ * softmax_alpha_t(m);
      for (int k = 0; k < alpha.size(); ++k) {
        if (m == k) {
          log_softmax_alpha(k).d_
              += alpha_ref.coeff(m).d_
                 + negative_alpha_m_d_times_softmax_alpha_t_m;
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
