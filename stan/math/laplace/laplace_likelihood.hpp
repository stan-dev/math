#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP

#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
namespace math {

// TO DO: create a parent structure, with each likelihood
// function acting as a child structure.
struct diff_poisson_log {
  Eigen::VectorXd n_samples_;
  Eigen::VectorXd sums_;

  diff_poisson_log();

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  template <typename T>
  T log_likelihood (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta)
    const {
    double factorial_term = 0;
    for (int i = 0; i < sums_.size(); i++)
      factorial_term += lgamma(sums_(i) + 1);

    return - factorial_term
      + dot_product(theta, sums_)
      - dot_product(exp(theta), n_samples_);
  }

  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    hessian = - elt_multiply(n_samples_, exp(theta));
    gradient = sums_ + hessian;
  }

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    return - elt_multiply(n_samples_, exp(theta));
  }
};

struct diff_logistic_log {
  Eigen::VectorXd n_samples_;
  Eigen::VectorXd sums_;

  diff_logistic_log();

  diff_logistic_log(const Eigen::VectorXd& n_samples,
                    const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  template <typename T>
  T log_likelihood (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta)
    const {
      Eigen::VectorXd one = rep_vector(1, theta.size());
      return sum(theta.cwiseProduct(sums_)
                   - n_samples_.cwiseProduct(log(one + exp(theta))));
    }

  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1> exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());

    gradient = sums_ - n_samples_.cwiseProduct(inv_logit(theta));

    hessian = - elt_multiply(n_samples_, elt_divide(exp_theta,
                                                    square(one + exp_theta)));
  }

  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    Eigen::VectorXd exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());
    Eigen::VectorXd nominator = exp_theta.cwiseProduct(exp_theta - one);
    Eigen::VectorXd denominator = square(one + exp_theta)
                                   .cwiseProduct(one + exp_theta);

    return n_samples_.cwiseProduct(elt_divide(nominator, denominator));
  }
};

// To experiment with the prototype, provide a built-in covariance
// function. In the final version, the user will pass the covariance
// function.
struct sqr_exp_kernel_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
           const T2& x, int M = 0) const {
    return stan::math::gp_exp_quad_cov(x, phi(0), phi(1))
    + 1e-9 * Eigen::MatrixXd::Identity(x.size(), x.size());
  }
};

}  // namespace math
}  // namespace stan

#endif
