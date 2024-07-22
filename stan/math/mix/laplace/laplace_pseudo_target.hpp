#ifndef STAN_MATH_LAPLACE_LAPLACE_PSEUDO_TARGET_HPP
#define STAN_MATH_LAPLACE_LAPLACE_PSEUDO_TARGET_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form.hpp>
#include <iostream>

namespace stan {
namespace math {
/**
 * Function to compute the pseudo target, $\tilde Z$,
 * with a custom derivative method.
 * NOTE: we actually don't need to compute the pseudo-target, only its
 * derivative.
 */
template <typename KMat,
          require_eigen_matrix_dynamic_vt<std::is_arithmetic, KMat>* = nullptr>
inline constexpr double laplace_pseudo_target(const KMat& K,
                                              const Eigen::VectorXd& a,
                                              const Eigen::MatrixXd& R,
                                              const Eigen::VectorXd& l_grad,
                                              const Eigen::VectorXd& s2) {
  // double s1 = 0.5 * quad_form(K, a) - 0.5 * sum((R * K).diagonal());
  // Eigen::VectorXd b = K * l_grad;
  // Eigen::VectorXd s3 = b - K * (R * b);
  // return s1 + s2.dot(s3);
  return 0;
}

/**
 * Overload function for case where K is passed as a matrix of var.
 */
template <typename KMat, typename AVec, typename RMat, typename LGradVec,
          typename S2Vec,
          require_eigen_matrix_dynamic_vt<is_var, KMat>* = nullptr>
inline auto laplace_pseudo_target(KMat&& K, AVec&& a, RMat&& R,
                                  LGradVec&& l_grad, S2Vec&& s2) {
  const Eigen::Index dim_theta = K.rows();
  auto K_arena = to_arena(std::forward<KMat>(K));
  auto&& a_ref = to_ref(std::forward<AVec>(a));
  auto&& R_ref = to_ref(std::forward<RMat>(R));
  auto&& s2_ref = to_ref(std::forward<S2Vec>(s2));
  auto&& l_grad_ref = to_ref(std::forward<LGradVec>(l_grad));
  arena_matrix<Eigen::MatrixXd> K_adj_arena
      = 0.5 * a_ref * a_ref.transpose() - 0.5 * R_ref
        + s2_ref * l_grad_ref.transpose()
        - (R_ref * (value_of(K_arena) * s2_ref)) * l_grad_ref.transpose();
  return make_callback_var(0.0, [K_arena, K_adj_arena](auto&& vi) mutable {
    K_arena.adj().array() += vi.adj() * K_adj_arena.array();
  });
}

}  // namespace math
}  // namespace stan

#endif
