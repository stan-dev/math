#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_NEG_BINOMIAL_2_LOG_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>

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
                          const std::vector<int>& y_index, int n_theta)
      : y_(y), y_index_(y_index), n_theta_(n_theta) {
    sums_ = Eigen::VectorXd::Zero(n_theta);
    n_samples_ = Eigen::VectorXd::Zero(n_theta);

    for (Eigen::Index i = 0; i < n_theta; i++) {
      n_samples_(y_index[i]) += 1;
      sums_(y_index[i]) += y[i];
    }
  }

  template <typename T_theta, typename T_eta>
  inline return_type_t<T_theta, T_eta> log_likelihood(
      const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    T_eta eta_scalar = eta(0);
    return_type_t<T_theta, T_eta> logp = 0;
    for (size_t i = 0; i < y_.size(); i++) {
      logp += binomial_coefficient_log(y_(i) + eta_scalar - 1, y_(i));
    }
    // CHECK -- is it better to vectorize this loop?
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    for (Eigen::Index i = 0; i < n_theta_; i++) {
      return_type_t<T_theta, T_eta> log_eta_plus_exp_theta
          = log(eta_scalar + exp_theta(i));
      logp += sums_(i) * (theta(i) - log_eta_plus_exp_theta)
              + n_samples_(i) * eta_scalar
                    * (log(eta_scalar) - log_eta_plus_exp_theta);
    }
    return logp;
  }

  template <typename T_theta, typename T_eta>
  inline Eigen::SparseMatrix<double> diff(
      const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
      Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>& gradient,
      const Eigen::Index hessian_block_size = 1) const {
    using scalar_t = return_type_t<T_theta, T_eta>;
    Eigen::VectorXd one = rep_vector(1, theta.size());
    const Eigen::Index theta_size = theta.size();
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1> sums_plus_n_eta
        = sums_ + eta_scalar * n_samples_;
    Eigen::Matrix<T_theta, Eigen::Dynamic, -1> exp_neg_theta = exp(-theta);

    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> one_plus_exp
        = one + eta_scalar * exp_neg_theta;
    gradient = sums_ - elt_divide(sums_plus_n_eta, one_plus_exp);
    Eigen::MatrixXd hessian_val = eta_scalar
                                  * sums_plus_n_eta.cwiseProduct(elt_divide(
                                      exp_neg_theta, square(one_plus_exp)));
    Eigen::SparseMatrix<double> hessian(theta_size, theta_size);
    hessian.reserve(Eigen::VectorXi::Constant(theta_size, hessian_block_size));
    // hessian.col(0) = - common_term;
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian.insert(i, i) = -hessian_val(i);
    }
    /*
        hessian = -eta_scalar
                  * sums_plus_n_eta.cwiseProduct(
                        elt_divide(exp_neg_theta, square(one_plus_exp)));
    */
    return hessian;
  }

  template <typename T_theta, typename T_eta>
  inline Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    using scalar_t = return_type_t<T_theta, T_eta>;
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1> eta_vec
        = rep_vector(eta_scalar, theta.size());
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> eta_plus_exp_theta
        = eta_vec + exp_theta;

    return -((sums_ + eta_scalar * n_samples_) * eta_scalar)
                .cwiseProduct(exp_theta.cwiseProduct(
                    elt_divide(eta_vec - exp_theta,
                               square(eta_plus_exp_theta)
                                   .cwiseProduct(eta_plus_exp_theta))));
  }

  template <typename T_theta, typename T_eta>
  inline Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    using scalar_t = return_type_t<T_theta, T_eta>;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_eta, Eigen::Dynamic, 1> y_plus_eta
        = y_ + rep_vector(eta_scalar, y_.size());
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> exp_theta_plus_eta
        = exp_theta + rep_vector(eta_scalar, theta.size());

    T_eta y_plus_eta_digamma_sum = 0;
    for (Eigen::Index i = 0; i < y_.size(); i++)
      y_plus_eta_digamma_sum += digamma(y_plus_eta(i));

    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> gradient_eta(1);
    gradient_eta(0)
        = y_plus_eta_digamma_sum - y_.size() * digamma(eta_scalar)
          - sum(elt_divide(sums_ + n_samples_ * eta_scalar, exp_theta_plus_eta))
          + sum(n_samples_ * log(eta_scalar)
                - n_samples_.cwiseProduct(log(exp_theta_plus_eta))
                + n_samples_);
    return gradient_eta;
  }

  template <typename T_theta, typename T_eta>
  inline auto diff_theta_eta(
      const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    using scalar_t = return_type_t<T_theta, T_eta>;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> diff_matrix(
        theta.size(), 1);
    diff_matrix.col(0) = -elt_divide(
        n_samples_ - sums_.cwiseProduct(exp_neg_theta),
        square(eta_scalar * exp_neg_theta + rep_vector(1, theta.size())));
    return diff_matrix;
  }

  // TODO: Address special case where we have an empty group (induces zero
  // elements in W).
  template <typename T_theta, typename T_eta>
  inline Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic,
                       Eigen::Dynamic>
  diff2_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) const {
    using scalar_t = return_type_t<T_theta, T_eta>;
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> one_plus_eta_exp
        = rep_vector(1, theta.size()) + eta_scalar * exp_neg_theta;

    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> diff_matrix(
        theta.size(), 1);

    diff_matrix.col(0) = -elt_divide(
        exp_neg_theta.cwiseProduct(-eta_scalar
                                       * exp_neg_theta.cwiseProduct(sums_)
                                   + sums_ + 2 * eta_scalar * n_samples_),
        square(one_plus_eta_exp).cwiseProduct(one_plus_eta_exp));  // );

    return diff_matrix;
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a negative binomial likelihood. Uses the 2nd parameterization.
 * Returns the marginal density p(y | phi) by marginalizing
 * out the latent gaussian variable, with a Laplace approximation.
 * See the laplace_marginal function for more details.
 *
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] y observations.
 * @param[in] y_index group to which each observation belongs. Each group
 *            is parameterized by one element of theta.
 * @param[in] covariance a function which returns the prior covariance.
 * @param[in] phi model parameters for the covariance functor.
 * @param[in] eta non-marginalized model parameters for the likelihood.
 * @param[in] x data for the covariance functor.
 * @param[in] delta additional real data for the covariance functor.
 * @param[in] delta_int additional integer data for covariance functor.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename CovarFun, typename Eta, typename Theta0, typename... Args>
inline auto laplace_marginal_tol_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    double tolerance, long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search, const Theta0& theta_0,
    CovarFun&& covariance_function, std::ostream* msgs, Args&&... args) {
  laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
  return laplace_marginal_density(
      diff_neg_binomial_2_log(to_vector(y), y_index, theta_0.size()),
      std::forward<CovarFun>(covariance_function), eta, theta_0, msgs,
      tolerance, max_num_steps, hessian_block_size, solver,
      max_steps_line_search, std::forward<Args>(args)...);
}

template <typename CovarFun, typename Eta, typename Theta0, typename... Args>
inline auto laplace_marginal_neg_binomial_2_log_lpmf(
    const std::vector<int>& y, const std::vector<int>& y_index, const Eta& eta,
    const Theta0& theta_0, CovarFun&& covariance_function, std::ostream* msgs,
    Args&&... args) {
  laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      diff_neg_binomial_2_log(to_vector(y), y_index, theta_0.size()),
      std::forward<CovarFun>(covariance_function), eta, theta_0, msgs,
      ops, std::forward<Args>(args)...);
}
}  // namespace math
}  // namespace stan

#endif
