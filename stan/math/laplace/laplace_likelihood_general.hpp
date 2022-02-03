#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP

// #include <stan/math/laplace/hessian_times_vector.hpp>
#include <stan/math/laplace/hessian_block_diag.hpp>

#include <Eigen/Sparse>

namespace stan {
namespace math {

/**
 * A structure to compute the log density, first, second,
 * and third-order derivatives for a likelihoood specified by the user.
 */
template <typename F>
struct diff_likelihood {
  /* Likelihood function. */
  F f_;
  /* Real variables passed to the likelihood. */
  Eigen::VectorXd delta_;
  /* Integer variables passed to the likelihood. */
  std::vector<int> delta_int_;
  /* stream to return print statements when function is called. */
  std::ostream* pstream_;
  template <typename FF, typename DeltaVec, typename DeltaInt>
  diff_likelihood(FF&& f, DeltaVec&& delta, DeltaInt&& delta_int,
                  std::ostream* pstream = 0)
      : f_(std::forward<FF>(f)),
        delta_(std::forward<DeltaVec>(delta)),
        delta_int_(std::forward<DeltaInt>(delta_int)),
        pstream_(pstream) {}

  template <typename T1, typename T2>
  inline T1 log_likelihood(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta) const {
    return f_(theta, eta, delta_, delta_int_, pstream_);
  }

  inline Eigen::SparseMatrix<double> diff(const Eigen::VectorXd& theta, const Eigen::VectorXd& eta,
                   Eigen::VectorXd& gradient,
                   const Eigen::Index hessian_block_size = 1) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;

    const Eigen::Index theta_size = theta.size();
    const Eigen::Index eta_size = eta.size();
    // CHECK -- do we need this scope?
    {
      nested_rev_autodiff nested;
      Matrix<var, Dynamic, 1> theta_var = theta;
      Matrix<var, Dynamic, 1> eta_var = eta;

      var f_var = f_(theta_var, eta_var, delta_, delta_int_, pstream_);
      grad(f_var.vi_);
      gradient.resize(theta_size + eta_size);
      for (Eigen::Index i = 0; i < theta_size; i++)
        gradient(i) = theta_var(i).adj();
      for (Eigen::Index i = 0; i < eta_size; i++)
        gradient(theta_size + i) = eta_var(i).adj();
    }

    if (hessian_block_size == 1) {
      Eigen::VectorXd v(theta_size);
      for (Eigen::Index i = 0; i < theta_size; i++)
        v(i) = 1;
      Eigen::VectorXd hessian_v = hessian_times_vector(f_, theta, eta, delta_,
                                                       delta_int_, v, pstream_);
      Eigen::SparseMatrix<double> hessian_theta(theta_size, theta_size);
      hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
      for (Eigen::Index i = 0; i < theta_size; i++)
        hessian_theta.insert(i, i) = hessian_v(i);
      return hessian_theta;
    } else {
      return hessian_block_diag(f_, theta, eta, delta_, delta_int_, hessian_block_size,
                         pstream_);
    }
  }

  inline Eigen::VectorXd third_diff(const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& eta) const {
    const Eigen::Index theta_size = theta.size();
    Eigen::VectorXd v = Eigen::VectorXd::Ones(theta_size);
    double f_theta;
    nested_rev_autodiff nested;
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta_var = theta;
    Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> theta_fvar(theta_size);
    for (Eigen::Index i = 0; i < theta_size; ++i) {
      theta_fvar(i) = fvar<var>(theta_var(i), v(i));
    }
    // TODO:(Steve) I think this can just be commented out?
    // fvar<var> ftheta_fvar = f_(theta_fvar, eta, delta_, delta_int_,
    // pstream_);

    Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> theta_ffvar(theta_size);
    for (Eigen::Index i = 0; i < theta_size; ++i) {
      theta_ffvar(i) = fvar<fvar<var>>(theta_fvar(i), v(i));
    }
    fvar<fvar<var>> ftheta_ffvar
        = f_(theta_ffvar, eta, delta_, delta_int_, pstream_);
    grad(ftheta_ffvar.d_.d_.vi_);
    return theta_var.adj();
  }

  inline Eigen::VectorXd compute_s2(const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& eta,
                                    const Eigen::MatrixXd& A,
                                    int hessian_block_size) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    nested_rev_autodiff nested;
    const Eigen::Index theta_size = theta.size();
    const Eigen::Index eta_size = eta.size();
    const Eigen::Index parm_size = theta_size + eta_size;
    // Matrix<var, Dynamic, 1> parm_var(parm_size);
    // for (Eigen::Index i = 0; i < theta_size; i++) parm_var(i) = theta(i);
    // for (Eigen::Index i = 0; i < eta_size; i++) parm_var(i + theta_size) =
    // eta(i);
    Matrix<var, Dynamic, 1> theta_var = theta;
    Matrix<var, Dynamic, 1> eta_var = eta;
    int n_blocks = theta_size / hessian_block_size;

    fvar<fvar<var>> target_ffvar = 0;

    for (Eigen::Index i = 0; i < hessian_block_size; ++i) {
      VectorXd v = VectorXd::Zero(theta_size);
      for (int j = i; j < theta_size; j += hessian_block_size)
        v(j) = 1;

      Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
      for (int j = 0; j < theta_size; ++j)
        theta_fvar(j) = fvar<var>(theta_var(j), v(j));

      Matrix<fvar<var>, Dynamic, 1> eta_fvar(eta_size);
      for (int j = 0; j < eta_size; ++j)
        eta_fvar(j) = fvar<var>(eta_var(j), 0);

      fvar<var> f_fvar = f_(theta_fvar, eta_fvar, delta_, delta_int_, pstream_);

      VectorXd w(theta_size);
      for (int j = 0; j < n_blocks; ++j) {
        for (int k = 0; k < hessian_block_size; ++k) {
          w(k + j * hessian_block_size)
              = A(k + j * hessian_block_size, i + j * hessian_block_size);
        }
      }

      Matrix<fvar<fvar<var>>, Dynamic, 1> theta_ffvar(theta_size);
      for (int j = 0; j < theta_size; ++j)
        theta_ffvar(j) = fvar<fvar<var>>(theta_fvar(j), w(j));

      Matrix<fvar<fvar<var>>, Dynamic, 1> eta_ffvar(eta_size);
      for (int j = 0; j < eta_size; ++j)
        eta_ffvar(j) = fvar<fvar<var>>(eta_fvar(j), 0);

      target_ffvar += f_(theta_ffvar, eta_ffvar, delta_, delta_int_, pstream_);
    }
    grad(target_ffvar.d_.d_.vi_);

    VectorXd parm_adj(parm_size);
    for (Eigen::Index i = 0; i < theta_size; ++i)
      parm_adj(i) = theta_var(i).adj();
    for (Eigen::Index i = 0; i < eta_size; ++i)
      parm_adj(theta_size + i) = eta_var(i).adj();

    return 0.5 * parm_adj;
  }

  inline Eigen::VectorXd diff_eta_implicit(const Eigen::VectorXd& v,
                                           const Eigen::VectorXd& theta,
                                           const Eigen::VectorXd& eta) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::VectorXd;

    nested_rev_autodiff nested;
    const Eigen::Index eta_size = eta.size();
    Matrix<var, Dynamic, 1> eta_var = eta;

    // CHECK -- can we avoid declaring theta as fvar<var>?
    // We currently compute derivatives wrt eta, which is not needed.
    const Eigen::Index theta_size = theta.size();
    Matrix<var, Dynamic, 1> theta_var = theta;
    Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
    for (Eigen::Index i = 0; i < theta_size; i++)
      theta_fvar(i) = fvar<var>(theta_var(i), v(i));

    // CHECK -- After merging develop branch, needed to do this.
    Matrix<fvar<var>, Dynamic, 1> eta_fvar(eta_size);
    for (Eigen::Index i = 0; i < eta_size; i++)
      eta_fvar(i) = fvar<var>(eta_var(i), 0);

    fvar<var> f_fvar = f_(theta_fvar, eta_fvar, delta_, delta_int_, pstream_);
    grad(f_fvar.d_.vi_);

    VectorXd diff_eta(eta_size);
    for (Eigen::Index i = 0; i < eta_size; i++)
      diff_eta(i) = eta_var(i).adj();
    return diff_eta;
  }
};

}  // namespace math
}  // namespace stan

#endif
