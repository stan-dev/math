#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP

// #include <stan/math/laplace/hessian_times_vector.hpp>
#include <stan/math/laplace/hessian_block_diag.hpp>
#include <stan/math/laplace/third_diff_directional.hpp>
#include <stan/math/laplace/partial_diff_theta.hpp>

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

  diff_likelihood(const F& f,
                  const Eigen::VectorXd& delta,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream = 0)
    : f_(f), delta_(delta), delta_int_(delta_int), pstream_(pstream) { }

  template <typename T1, typename T2>
  T1 log_likelihood(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
                    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta)
    const {
      return f_(theta, eta, delta_, delta_int_, pstream_);
    }

  void diff (const Eigen::VectorXd& theta,
             const Eigen::VectorXd& eta,
             Eigen::VectorXd& gradient,
             Eigen::SparseMatrix<double>& hessian_theta,
             int hessian_block_size = 1) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;

    int theta_size = theta.size();
    int eta_size = eta.size();
    // CHECK -- do we need this scope?
    {
      nested_rev_autodiff nested;
      Matrix<var, Dynamic, 1> theta_var = theta;
      Matrix<var, Dynamic, 1> eta_var = eta;

      var f_var = f_(theta_var, eta_var, delta_, delta_int_, pstream_);
      grad(f_var.vi_);
      gradient.resize(theta_size + eta_size);
      for (int i = 0; i < theta_size; i++) gradient(i) = theta_var(i).adj();
      for (int i = 0; i < eta_size; i++)
        gradient(theta_size + i) = eta_var(i).adj();
    }

    hessian_theta.resize(theta_size, theta_size);
    double f_theta;
    if (hessian_block_size == 1) {
      Eigen::VectorXd v(theta_size);
      for (int i = 0; i < theta_size; i++) v(i) = 1;
      Eigen::VectorXd hessian_v;
      hessian_times_vector(f_, theta, eta, delta_, delta_int_,
                           v, f_theta, hessian_v, pstream_);
      hessian_theta.reserve(Eigen::VectorXi::Constant(theta_size, 1));
      for (int i = 0; i < theta_size; i++)
        hessian_theta.insert(i, i) = hessian_v(i);
    } else {
      hessian_block_diag(f_, theta, eta, delta_, delta_int_,
                         hessian_block_size,
                         f_theta, hessian_theta, pstream_);
    }
  }

  Eigen::VectorXd third_diff(const Eigen::VectorXd& theta,
                             const Eigen::VectorXd& eta) const {

    int theta_size = theta.size();
    Eigen::VectorXd v(theta_size);
    for (int i = 0; i < theta_size; i++) v(i) = 1;
    double f_theta;
    Eigen::VectorXd third_diff_tensor;

    third_diff_directional(f_, theta, eta, delta_, delta_int_,
                           f_theta, third_diff_tensor,
                           v, v, pstream_);

    return third_diff_tensor;
  }

  Eigen::VectorXd compute_s2(const Eigen::VectorXd& theta,
                             const Eigen::VectorXd& eta,
                             const Eigen::MatrixXd& A,
                             int hessian_block_size) const {
    return partial_diff_theta(f_, theta, eta, delta_, delta_int_, A,
                              hessian_block_size, pstream_);
  }

  Eigen::VectorXd diff_eta_implicit(const Eigen::VectorXd& v,
                                    const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& eta) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using Eigen::VectorXd;

    nested_rev_autodiff nested;
    int eta_size = eta.size();
    Matrix<var, Dynamic, 1> eta_var = eta;

    // CHECK -- can we avoid declaring theta as fvar<var>?
    // We currently compute derivatives wrt eta, which is not needed.
    int theta_size = theta.size();
    Matrix<var, Dynamic, 1> theta_var = theta;
    Matrix<fvar<var>, Dynamic, 1> theta_fvar(theta_size);
    for (int i = 0; i < theta_size; i++)
      theta_fvar(i) = fvar<var>(theta_var(i), v(i));

    // CHECK -- After merging develop branch, needed to do this.
    Matrix<fvar<var>, Dynamic, 1> eta_fvar(eta_size);
    for (int i = 0; i < eta_size; i++) eta_fvar(i) = fvar<var>(eta_var(i), 0);

    fvar<var> f_fvar = f_(theta_fvar, eta_fvar, delta_, delta_int_, pstream_);
    grad(f_fvar.d_.vi_);

    VectorXd diff_eta(eta_size);
    for (int i = 0; i < eta_size; i++) diff_eta(i) = eta_var(i).adj();
    return diff_eta;
  }


  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    std::cout << "THIS FUNCTION SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1> void_matrix;
    return void_matrix;
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    std::cout << "THIS FUNCTION SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic,
      Eigen::Dynamic> void_matrix;
    return void_matrix;
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff2_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta)
  const {
    std::cout << "THIS FUNCTION SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic,
      Eigen::Dynamic> void_matrix;
    return void_matrix;
  }
};

}  // namespace math
}  // namespace stan

#endif
