#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_POISSON_LOG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_POISSON_LOG_HPP

#include <stan/math/prim/fun/lgamma.hpp>

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
   * @param[in] eta_dummy additional parameters (use for other likelihoods).
   * @return the log density.
   */
  template <typename T1, typename T2>
  T1 log_likelihood (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
                     const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy)
    const {
    double factorial_term = 0;
    for (int i = 0; i < sums_.size(); i++)
      factorial_term += lgamma(sums_(i) + 1);
    Eigen::Matrix<T1, Eigen::Dynamic, 1> shifted_mean = theta + log_exposure_;

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
   * @param[in] eta_dummy additional parameters (use for other likelihoods).
   * @param[in, out] gradient
   * @param[in, out] hessian diagonal, so stored in a vector.
   */
  template <typename T1, typename T2>
  void diff (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy,
             Eigen::Matrix<T1, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T1, Eigen::Dynamic, 1>& hessian) const {
    Eigen::Matrix<T1, Eigen::Dynamic, 1>
      common_term = n_samples_.cwiseProduct(exp(theta + log_exposure_));

    gradient = sums_ - common_term;
    hessian = - common_term;
  }

  /**
   * Returns the third derivative tensor. Because it is ("cubic") diagonal,
   * the object is stored in a vector.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @param[in] eta_dummy additional parameters (use for other likelihoods).
   * @return A vector containing the non-zero elements of the third
   *         derivative tensor.
   */
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy) const {
    return -n_samples_.cwiseProduct(exp(theta + log_exposure_));
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
