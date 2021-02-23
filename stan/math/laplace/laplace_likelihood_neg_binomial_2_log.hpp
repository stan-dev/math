#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_NEG_BINOMIAL_2_LOG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_NEG_BINOMIAL_2_LOG_HPP

#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct diff_neg_binomial_2_log {
  /* Observed counts */
  Eigen::VectorXd y_;
  /* Latent parameter index for each observation. */
  std::vector<int> y_index_;
  /* The number of samples in each group. */
  Eigen::VectorXd n_samples_;
  /* The sum of cours in each group. */
  Eigen::VectorXd sums_;
  /* Number of latent Gaussian variables. */
  int n_theta_;

  diff_neg_binomial_2_log(const Eigen::VectorXd& y,
                          const std::vector<int>& y_index,
                          int n_theta)
    : y_(y), y_index_(y_index), n_theta_(n_theta) {
    sums_ = Eigen::VectorXd::Zero(n_theta);
    n_samples_ = Eigen::VectorXd::Zero(n_theta);

    for (int i = 0; i < n_theta; i++) {
      n_samples_(y_index[i]) += 1;
      sums_(y_index[i]) += y[i];
    }
  }

  template <typename T_theta, typename T_eta>
  return_type_t<T_theta, T_eta>
  log_likelihood (const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    T_eta eta_scalar = eta(0);
    return_type_t<T_theta, T_eta> logp = 0;
    for (size_t i = 0; i < y_.size(); i++) {
      logp += binomial_coefficient_log(y_(i) + eta_scalar - 1, y_(i));
    }
    // CHECK -- is it better to vectorize this loop?
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    for (int i = 0; i < n_theta_; i++) {
      return_type_t<T_theta, T_eta>
        log_eta_plus_exp_theta = log(eta_scalar + exp_theta(i));
      logp += sums_(i) * (theta(i) - log_eta_plus_exp_theta)
               + n_samples_(i) * eta_scalar
               * (log(eta_scalar) - log_eta_plus_exp_theta);
    }
    return logp;
  }

  template <typename T_theta, typename T_eta>
  void diff (const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
             Eigen::Matrix<return_type_t<T_theta, T_eta>,
               Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<return_type_t<T_theta, T_eta>,
               Eigen::Dynamic, 1>& hessian) const {
    typedef return_type_t<T_theta, T_eta> scalar;
    Eigen::VectorXd one = rep_vector(1, theta.size());
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1>
      sums_plus_n_eta = sums_ + eta_scalar * n_samples_;
    Eigen::Matrix<T_theta, Eigen::Dynamic, -1> exp_neg_theta = exp(-theta);

    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      one_plus_exp = one + eta_scalar * exp_neg_theta;
    gradient = sums_ - elt_divide(sums_plus_n_eta, one_plus_exp);

    hessian = - eta_scalar * sums_plus_n_eta.
      cwiseProduct(elt_divide(exp_neg_theta, square(one_plus_exp)));
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    typedef return_type_t<T_theta, T_eta> scalar;
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1>
      eta_vec = rep_vector(eta_scalar, theta.size());
    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      eta_plus_exp_theta = eta_vec + exp_theta;

    return - ((sums_ + eta_scalar * n_samples_) * eta_scalar).
      cwiseProduct(exp_theta.cwiseProduct(
        elt_divide(eta_vec - exp_theta,
        square(eta_plus_exp_theta).cwiseProduct(eta_plus_exp_theta))));
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    typedef return_type_t<T_theta, T_eta> scalar;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1>
      y_plus_eta = y_ + rep_vector(eta_scalar, y_.size());
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    Eigen::Matrix<scalar, Eigen::Dynamic, 1>
      exp_theta_plus_eta = exp_theta + rep_vector(eta_scalar, theta.size());

    T_eta y_plus_eta_digamma_sum = 0;
    for (int i = 0; i < y_.size(); i++)
      y_plus_eta_digamma_sum += digamma(y_plus_eta(i));

     Eigen::Matrix<scalar, Eigen::Dynamic, 1> gradient_eta(1);
     gradient_eta(0) =
       y_plus_eta_digamma_sum - y_.size() * digamma(eta_scalar)
        - sum(elt_divide(sums_ + n_samples_ * eta_scalar, exp_theta_plus_eta))
        + sum(n_samples_ * log(eta_scalar)
              - n_samples_.cwiseProduct(log(exp_theta_plus_eta))
              + n_samples_);
      return gradient_eta;
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    typedef return_type_t<T_theta, T_eta> scalar;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);
    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
      diff_matrix(theta.size(), 1);
    diff_matrix.col(0)
      = - elt_divide(n_samples_ - sums_.cwiseProduct(exp_neg_theta),
      square(eta_scalar * exp_neg_theta + rep_vector(1, theta.size())));
    return diff_matrix;
  }

  // TODO: Address special case where we have an empty group (induces zero
  // elements in W).
  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, Eigen::Dynamic>
  diff2_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta)
  const {
    typedef return_type_t<T_theta, T_eta> scalar;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> one_plus_eta_exp
      = rep_vector(1, theta.size()) + eta_scalar * exp_neg_theta;

    Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic>
      diff_matrix(theta.size(), 1);

    diff_matrix.col(0) =
      - elt_divide(exp_neg_theta.cwiseProduct(
      - eta_scalar * exp_neg_theta.cwiseProduct(sums_)
      + sums_ + 2 * eta_scalar * n_samples_),
      square(one_plus_eta_exp).cwiseProduct(one_plus_eta_exp)); // );

    return diff_matrix;
  }
};

}  // namespace math
}  // namespace stan

#endif
