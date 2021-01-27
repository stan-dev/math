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
      sums_(y_index[i]) += y[i];
    }
  }

  template <typename T_theta, typename T_eta>
  return_type_t<T_theta, T_eta>
  log_likelihood (const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) {
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
             const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) {
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
           const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) {
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
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                 const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta) {
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);

    return - elt_divide(n_samples_ - sums_.cwiseProduct(exp_neg_theta),
      square(eta_scalar * exp_neg_theta + rep_vector(1, theta.size())));
  }

  template <typename T_theta, typename T_eta>
  Eigen::Matrix<return_type_t<T_theta, T_eta>, Eigen::Dynamic, 1>
  diff2_theta_eta(const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& theta,
                  const Eigen::Matrix<T_eta, Eigen::Dynamic, 1>& eta,
                  const Eigen::Matrix<T_theta, Eigen::Dynamic, 1>& W_root) {
    T_eta eta_scalar = eta(0);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> exp_neg_theta = exp(-theta);
    Eigen::Matrix<T_theta, Eigen::Dynamic, 1> one_plus_eta_exp
      = rep_vector(1, theta.size()) + eta_scalar * exp_neg_theta;

    return 0.5 * (W_root.cwiseInverse()).cwiseProduct(
      elt_divide(exp_neg_theta.cwiseProduct(
      - eta_scalar * exp_neg_theta.cwiseProduct(sums_)
      + sums_ + 2 * eta_scalar * n_samples_),
      square(one_plus_eta_exp).cwiseProduct(one_plus_eta_exp)));
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
