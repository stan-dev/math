#ifndef STAN_MATH_LAPLACE_PARTIAL_DIFF_THETA_HPP
#define STAN_MATH_LAPLACE_PARTIAL_DIFF_THETA_HPP

// TODO: refine include.
#include <stan/math/prim/fun/col.hpp>
#include <stan/math/mix.hpp>

namespace stan {
namespace math {
  /**
   * Returns the partial derivative of the approximate marginal
   * distribution with respect to theta and eta.
   * The derivative with respect to theta is denoted s2 in
   * laplace_marginal.hpp.
   */
   // TODO: rename function, since we also differentiate wrt eta.
   // TODO: address case where eta / theta are doubles and we don't
   // want full derivatives.
   template <typename F>
   Eigen::VectorXd partial_diff_theta(const F& f,
                                      const Eigen::VectorXd& theta,
                                      const Eigen::VectorXd& eta,
                                      const Eigen::VectorXd& delta,
                                      const std::vector<int>& delta_int,
                                      const Eigen::MatrixXd& A,
                                      int hessian_block_size,
                                      std::ostream* pstream = 0) {
    using Eigen::VectorXd;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::Dynamic;

    nested_rev_autodiff nested;
    int theta_size = theta.size();
    int eta_size = eta.size();
    int parm_size = theta_size + eta_size;
    // Matrix<var, Dynamic, 1> parm_var(parm_size);
    // for (int i = 0; i < theta_size; i++) parm_var(i) = theta(i);
    // for (int i = 0; i < eta_size; i++) parm_var(i + theta_size) = eta(i);
    Matrix<var, Dynamic, 1> theta_var = theta;
    Matrix<var, Dynamic, 1> eta_var = eta;
    int n_blocks = theta_size / hessian_block_size;

    fvar<fvar<var>> target_ffvar = 0;

    for (int i = 0; i < hessian_block_size; ++i) {
      VectorXd v = VectorXd::Zero(theta_size);
      for (int j = i; j < theta_size; j += hessian_block_size) v(j) = 1;

      Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
      for (int j = 0; j < theta_size; ++j)
        theta_fvar(j) = fvar<var>(theta_var(j), v(j));

      Matrix<fvar<var>, Dynamic, 1> eta_fvar(eta_size);
      for (int j = 0; j < eta_size; ++j) eta_fvar(j) = fvar<var>(eta_var(j), 0);

      fvar<var> f_fvar = f(theta_fvar, eta_fvar, delta, delta_int, pstream);

      VectorXd w(theta_size);
      for (int j = 0; j < n_blocks; ++j) {
        for (int k = 0; k < hessian_block_size; ++k) {
          w(k + j * hessian_block_size) =  A(k + j * hessian_block_size,
                                             i + j * hessian_block_size);
        }
      }

      Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
      for (int j = 0; j < theta_size; ++j)
        theta_ffvar(j) = fvar<fvar<var>>(theta_fvar(j), w(j));

      Matrix<fvar<fvar<var>>, Dynamic, 1> eta_ffvar(eta_size);
      for (int j = 0; j < eta_size; ++j)
        eta_ffvar(j) = fvar<fvar<var>>(eta_fvar(j), 0);

      target_ffvar +=
        f(theta_ffvar, eta_ffvar, delta, delta_int, pstream);
    }
    grad(target_ffvar.d_.d_.vi_);

    VectorXd parm_adj(parm_size);
    for (int i = 0; i < theta_size; ++i) parm_adj(i) = theta_var(i).adj();
    for (int i = 0; i < eta_size; ++i)
      parm_adj(theta_size + i) = eta_var(i).adj();

    return 0.5 * parm_adj;

    // VectorXd theta_adj(theta_size);
    // for (int i = 0; i < theta_size; ++i) theta_adj(i) = theta_var(i).adj();
    // return 0.5 * theta_adj;
  }

}  // namespace math
}  // namespace stan

#endif
