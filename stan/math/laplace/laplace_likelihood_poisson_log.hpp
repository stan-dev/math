#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_POISSON_LOG_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_POISSON_LOG_HPP

#include <stan/math/prim/fun/lgamma.hpp>
#include <Eigen/Sparse>

namespace stan {
namespace math {

struct poisson_log_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * @tparam T_theta Type of the log Poisson rate.
   * @tparam T_eta Type of the auxiliary parameter (not used here).
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y sum of counts in each group.
   * @param[in] delta_int number of observations in each group.
   * return lpmf for a Poisson with a log link.
   */
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, -1, 1>& theta,
      const Eigen::Matrix<T_eta, -1, 1>& eta, const Eigen::VectorXd& y,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    Eigen::VectorXd n_samples = to_vector(delta_int);
    return -lgamma(y.array() + 1).sum() + theta.dot(y)
           - n_samples.dot(exp(theta));
  }
};

struct poisson_log_exposure_likelihood {
  /**
   * Returns the lpmf for a Poisson with a log link across
   * multiple groups. No need to compute the log normalizing constant.
   * Same as above, but includes a exposure term to correct the
   * log rate for each group.
   * @tparam T_theta Type of the log Poisson rate.
   * @tparam T_eta Type of the auxiliary parameter (not used here).
   * @param[in] theta log Poisson rate for each group.
   * @param[in] y_and_ye First n elements contain the sum of counts
   *                     in each group, next n elements the exposure
   *                     in each group, where n is the number of groups.
   * @param[in] delta_int number of observations in each group.
   * return lpmf for a Poisson with a log link.
   */
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, -1, 1>& theta,
      const Eigen::Matrix<T_eta, -1, 1>& eta, const Eigen::VectorXd& y_and_ye,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    int n = delta_int.size();
    Eigen::VectorXd y = y_and_ye.head(n);
    Eigen::VectorXd ye = y_and_ye.tail(n);

    Eigen::VectorXd n_samples = to_vector(delta_int);
    Eigen::Matrix<T_theta, -1, 1> shifted_mean = theta + log(ye);
    return -lgamma(y.array() + 1).sum() + shifted_mean.dot(y)
           - n_samples.dot(exp(shifted_mean));
  }
};

// NOTE: might not need the code below anymore, since we'll switch to
// the general diff likelihood.
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

  template <typename SampleVec, typename SumVec,
            require_all_eigen_vector_t<SampleVec, SumVec>* = nullptr>
  diff_poisson_log(SampleVec&& n_samples, SumVec&& sums)
      : n_samples_(std::forward<SampleVec>(n_samples)),
        sums_(std::forward<SumVec>(sums)),
        log_exposure_(Eigen::VectorXd::Zero(sums_.size())) {}

  template <typename SampleVec, typename SumVec, typename LogExpVec,
            require_all_eigen_vector_t<SampleVec, SumVec, LogExpVec>* = nullptr>
  diff_poisson_log(SampleVec&& n_samples, SumVec&& sums,
                   LogExpVec&& log_exposure)
      : n_samples_(std::forward<SampleVec>(n_samples)),
        sums_(std::forward<SumVec>(sums)),
        log_exposure_(std::forward<LogExpVec>(log_exposure)) {}

  /**
   * Return the log density.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @param[in] eta_dummy additional parameters (use for other likelihoods).
   * @return the log density.
   */
  template <typename T1, typename T2>
  inline auto log_likelihood(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy) const {
    double factorial_term = 0;
    for (Eigen::Index i = 0; i < sums_.size(); i++)
      factorial_term += lgamma(sums_(i) + 1);
    Eigen::Matrix<T1, Eigen::Dynamic, 1> shifted_mean = theta + log_exposure_;

    // return shifted_mean.dot(sums_) - n_samples_.dot(exp(shifted_mean));
    return -lgamma(sums_.array() + 1).sum() + (shifted_mean).dot(sums_)
           - n_samples_.dot(exp(shifted_mean));
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
  inline Eigen::SparseMatrix<double> diff(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy,
      Eigen::Matrix<T1, Eigen::Dynamic, 1>& gradient,
      // Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& hessian,
      const Eigen::Index hessian_block_size = 1) const {
    const Eigen::Index theta_size = theta.size();
    Eigen::Matrix<T1, Eigen::Dynamic, 1> common_term
        = n_samples_.cwiseProduct(exp(theta + log_exposure_));

    gradient = sums_ - common_term;
    Eigen::SparseMatrix<double> hessian(theta_size, theta_size);
    hessian.reserve(Eigen::VectorXi::Constant(theta_size, hessian_block_size));
    // hessian.col(0) = - common_term;
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian.insert(i, i) = -common_term(i);
    }
    return hessian;
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
  inline Eigen::Matrix<T1, Eigen::Dynamic, 1> third_diff(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy) const {
    return -n_samples_.cwiseProduct(exp(theta + log_exposure_));
  }
};

}  // namespace math
}  // namespace stan

#endif
