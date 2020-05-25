#ifndef STAN_MATH_FWD_FUN_SOFTMAX_HPP
#define STAN_MATH_FWD_FUN_SOFTMAX_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/softmax.hpp>

namespace stan {
namespace math {

template <typename ColVec,
          require_eigen_col_vector_vt<is_fvar, ColVec>* = nullptr>
inline auto softmax(const ColVec& alpha) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T = typename value_type_t<ColVec>::Scalar;

  Matrix<T, Dynamic, 1> alpha_t(alpha.size());
  for (int k = 0; k < alpha.size(); ++k) {
    alpha_t.coeffRef(k) = alpha.coeff(k).val_;
  }

  Matrix<T, Dynamic, 1> softmax_alpha_t = softmax(alpha_t);

  Matrix<fvar<T>, Dynamic, 1> softmax_alpha(alpha.size());
  for (int k = 0; k < alpha.size(); ++k) {
    softmax_alpha.coeffRef(k).val_ = softmax_alpha_t.coeff(k);
    softmax_alpha.coeffRef(k).d_ = 0;
  }

  for (int m = 0; m < alpha.size(); ++m) {
    T negative_alpha_m_d_times_softmax_alpha_t_m
        = -alpha.coeff(m).d_ * softmax_alpha_t.coeff(m);
    for (int k = 0; k < alpha.size(); ++k) {
      if (m == k) {
        softmax_alpha.coeffRef(k).d_
            += softmax_alpha_t.coeff(k)
               * (alpha.coeff(m).d_
                  + negative_alpha_m_d_times_softmax_alpha_t_m);
      } else {
        softmax_alpha.coeffRef(k).d_
            += softmax_alpha_t.coeff(k)
               * negative_alpha_m_d_times_softmax_alpha_t_m;
      }
    }
  }

  return softmax_alpha;
}

}  // namespace math
}  // namespace stan
#endif
