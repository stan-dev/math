#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP

#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
namespace math {

// TO DO: create a parent structure, with each likelihood
// function acting as a child structure.

/**
 * Create an Eigen vector whose elements are all ones.
 */
Eigen::VectorXd init_one(int n) {
  Eigen::VectorXd ones(n);
  for (int i = 0; i < n; i++) ones(i) = 1;
  return ones;
}

/**
 * A structure to compute the log density, first, second,
 * and third-order derivatives for a log poisson likelihood
 * whith multiple groups.
 * This structure can be passed to the the laplace_marginal function.
 * Uses sufficient statistics for the data. 
 */
struct diff_poisson_log {
  /* The number of samples in each group. */
  Eigen::VectorXd n_samples_;
  /* The sum of counts in each group. */
  Eigen::VectorXd sums_;
  /* exposure, i.e. off-set term for the latent variable. */
  Eigen::VectorXd exposure_;

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) {
    exposure_ = init_one(sums_.size());
  }

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums,
                   const Eigen::VectorXd& exposure)
    : n_samples_(n_samples), sums_(sums), exposure_(exposure) { }
  
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

    return - factorial_term
      + (exposure_.cwiseProduct(theta)).dot(sums_) -
        exp(exposure_.cwiseProduct(theta)).dot(n_samples_);
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
      common_term = - n_samples_.cwiseProduct(exposure_)
                        .cwiseProduct(exp(theta.cwiseProduct(exposure_)));
    hessian = common_term.cwiseProduct(exposure_);
    gradient = sums_.cwiseProduct(exposure_) + common_term;
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
    return - n_samples_.cwiseProduct(exp(theta.cwiseProduct(exposure_)))
                       .cwiseProduct(exposure_)
                       .cwiseProduct(exposure_)
                       .cwiseProduct(exposure_);
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

// TO DO: delete this structure.
// To experiment with the prototype, provide a built-in covariance
// function. In the final version, the user will pass the covariance
// function.
struct sqr_exp_kernel_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
           const T2& x, int M = 0) const {
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
