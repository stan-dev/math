#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP

#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

// TO DO: create a parent structure, with each likelihood
// function acting as a child structure.

/**
 * A structure to compute the log density, first, second,
 * and third-order derivatives for a log poisson likelihood
 * whith multiple groups.
 * This structure can be passed to the the laplace_marginal function.
 * Uses sufficient statistics for the data.
 */
 // FIX ME -- cannot use the sufficient statistic to compute log density in
 // because of log factorial term.
struct diff_poisson_log {
  /* The number of samples in each group. */
  Eigen::VectorXd n_samples_;
  /* The sum of counts in each group. */
  Eigen::VectorXd sums_;
  /* exposure, i.e. off-set term for the latent variable. */
  Eigen::VectorXd log_exposure_;

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) {
    log_exposure_ = Eigen::VectorXd::Zero(sums.size());
  }

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums,
                   const Eigen::VectorXd& log_exposure)
    : n_samples_(n_samples), sums_(sums), log_exposure_(log_exposure) { }

  /**
   * Return the log density.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return the log density.
   */
  template <typename T>
  T log_likelihood (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta)
    const {
    double factorial_term = 0;
    for (int i = 0; i < sums_.size(); i++)
      factorial_term += lgamma(sums_(i) + 1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> shifted_mean = theta + log_exposure_;

    return - factorial_term
      + (shifted_mean).dot(sums_) - n_samples_.dot(exp(shifted_mean));
  }

  /**
   * Returns the gradient of the log density, and the hessian.
   * Since the latter is diagonal, it is stored inside a vector.
   * The two objects are computed together, because we always use
   * both when solving the Newton iteration of the Laplace
   * approximation, and to avoid redundant computation.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @param[in, out] gradient
   * @param[in, out] hessian diagonal, so stored in a vector.
   */
  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1>
      common_term = n_samples_.cwiseProduct(exp(theta + log_exposure_));

    gradient = sums_ - common_term;
    hessian = - common_term;
  }

  /**
   * Returns the third derivative tensor. Because it is ("cubic") diagonal,
   * the object is stored in a vector.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return A vector containing the non-zero elements of the third
   *         derivative tensor.
   */
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    return -n_samples_.cwiseProduct(exp(theta + log_exposure_));
  }
};

/**
 * A structure to compute the log density, first, second,
 * and third-order derivatives for a Bernoulli logistic likelihood
 * whith multiple groups.
 * This structure can be passed to the the laplace_marginal function.
 * Uses sufficient statistics for the data.
 */
struct diff_logistic_log {
  /* The number of samples in each group. */
  Eigen::VectorXd n_samples_;
  /* The sum of counts in each group. */
  Eigen::VectorXd sums_;

  diff_logistic_log(const Eigen::VectorXd& n_samples,
                    const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  /**
   * Return the log density.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return the log density.
   */
  template <typename T>
  T log_likelihood (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta)
    const {
      Eigen::VectorXd one = rep_vector(1, theta.size());
      return sum(theta.cwiseProduct(sums_)
                   - n_samples_.cwiseProduct(log(one + exp(theta))));
    }

  /**
   * Returns the gradient of the log density, and the hessian.
   * Since the latter is diagonal, it is stored inside a vector.
   * The two objects are computed together, because we always use
   * both when solving the Newton iteration of the Laplace
   * approximation, and to avoid redundant computation.
   * @tparam T type of the Bernoulli logistic parameter.
   * @param[in] theta Bernoulli logistic parameters for each group.
   * @param[in, out] gradient
   * @param[in, out] hessian diagonal, so stored in a vector.
   */
  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1> exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());

    gradient = sums_ - n_samples_.cwiseProduct(inv_logit(theta));

    hessian = - n_samples_.cwiseProduct(elt_divide(exp_theta,
                                                    square(one + exp_theta)));
  }

  /**
   * Returns the third derivative tensor. Because it is (cubic) diagonal,
   * the object is stored in a vector.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return A vector containing the non-zero elements of the third
   *         derivative tensor.
   */
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
      sums_(y_index[i]) += 1;
    }
  }

  template <typename T_theta, typename T_eta>
  return_type_t<T_theta, T_eta>
  log_likelihood (const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const T_eta& eta) {
    return_type_t<T_theta, T_eta> logp = 0;
    for (size_t i = 0; i < y_.size(); i++) {
      logp += binomial_coefficient_log(y_(i) + eta - 1, y_(i));
    }
    // CHECK -- is it better to vectorize this loop?
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_theta = exp(theta);
    for (int i = 0; i < n_theta_; i++) {
      return_type_t<T_theta, T_eta>
        log_theta_plus_exp_theta = log(exp_theta(i) + eta);
      logp += y_(i) * (theta(i) - log_theta_plus_exp_theta)
               + n_samples_(i) * eta * (log(eta) - log_theta_plus_exp_theta);
    }
    return logp;
  }
};

struct diff_student_t {
  /* Observations. */
  Eigen::VectorXd y_;
  /* Latent parameter index for each observation. */
  std::vector<int> y_index_;
  // QUESTION - Save eta here too?

  diff_student_t(const Eigen::VectorXd& y,
                 const std::vector<int>& y_index)
    : y_(y), y_index_(y_index) { }

  /**
   * Returns the log density.
   */
  template <typename T_theta, typename T_eta>
  return_type_t<T_theta, T_eta>
  log_likelihood (const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta)
    const {
    T_eta nu = eta(0);
    T_eta sigma = eta(1);
    T_eta sigma_squared = sigma * sigma;

    int n = theta.size();

    // CHECK -- probably don't need normalizing constant.
    return_type_t<T_theta, T_eta>
    log_constant = n * (lgamma((nu + 1) / 2) - lgamma(nu / 2)
                    - LOG_SQRT_PI - 0.5 * log(nu) - log(sigma));

    T_theta log_kernel = 0;

    for (int i = 0; i < n; i++) {
      T_theta distance = y_(i) - theta(y_index_[i]);
      log_kernel += log(1 + distance * distance / (nu * sigma_squared));
    }

    return log_constant - 0.5 * (nu + 1) * log_kernel;
  }
};


// TO DO: delete this structure.
// To experiment with the prototype, provide a built-in covariance
// function. In the final version, the user will pass the covariance
// function.
struct sqr_exp_kernel_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
              const T2& x,
              const std::vector<double>& delta,
              const std::vector<int>& delta_int,
              std::ostream* msgs = nullptr) const {
    double jitter = 1e-8;
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
      kernel = stan::math::gp_exp_quad_cov(x, phi(0), phi(1));
    for (int i = 0; i < kernel.cols(); i++)
      kernel(i, i) += jitter;

    return kernel;
  }
};



}  // namespace math
}  // namespace stan

#endif
