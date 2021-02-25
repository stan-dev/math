#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_GENERAL_HPP

#include <stan/math/laplace/hessian_times_vector.hpp>
#include <stan/math/laplace/third_diff_directional.hpp>

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
             Eigen::VectorXd& hessian) const {
    using Eigen::Matrix;
    using Eigen::Dynamic;

    int theta_size = theta.size();
    // CHECK -- do we need this scope?
    {
      nested_rev_autodiff nested;
      Matrix<var, Dynamic, 1> theta_var = theta;
      var f_var = f_(theta_var, eta, delta_, delta_int_, pstream_);
      grad(f_var.vi_);
      gradient.resize(theta_size);
      for (int i = 0; i < theta_size; i++) gradient(i) = theta_var(i).adj();
    }

    Eigen::VectorXd v(theta_size);
    for (int i = 0; i < theta_size; i++) v(i) = 1;
    double f_theta;
    hessian_times_vector(f_, theta, eta, delta_, delta_int_,
                         v, f_theta, hessian, pstream_);
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

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    std::cout << "THIS FUNCTIONS SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1> void_matrix;
    return void_matrix;
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    std::cout << "THIS FUNCTIONS SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic,
      Eigen::Dynamic> void_matrix;
    return void_matrix;
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff2_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta)
  const {
    std::cout << "THIS FUNCTIONS SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic,
      Eigen::Dynamic> void_matrix;
    return void_matrix;
  }
};

}  // namespace math
}  // namespace stan

#endif
