#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_BERNOULLI_LOGIT_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_BERNOULLI_LOGIT_HPP


namespace stan {
namespace math {

/**
 * A structure to compute the log density, first, second,
 * and third-order derivatives for a Bernoulli logistic likelihood
 * whith multiple groups.
 * This structure can be passed to the the laplace_marginal function.
 * Uses sufficient statistics for the data.
 */
struct diff_bernoulli_logit {
  /* The number of samples in each group. */
  Eigen::VectorXd n_samples_;
  /* The sum of counts in each group. */
  Eigen::VectorXd sums_;

  diff_bernoulli_logit(const Eigen::VectorXd& n_samples,
                       const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  /**
   * Return the log density.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return the log density.
   */
  template <typename T1, typename T2>
  T1 log_likelihood (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
                    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy)
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
  template <typename T1, typename T2>
  void diff (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy,
             Eigen::Matrix<T1, Eigen::Dynamic, 1>& gradient,
             Eigen::SparseMatrix<double>& hessian,
             // Eigen::Matrix<T1, Eigen::Dynamic, 1>& hessian,
             int block_size_dummy) const {
    Eigen::Matrix<T1, Eigen::Dynamic, 1> exp_theta = exp(theta);
    int theta_size = theta.size();
    Eigen::VectorXd one = rep_vector(1, theta_size);

    gradient = sums_ - n_samples_.cwiseProduct(inv_logit(theta));

    Eigen::Matrix<T1, Eigen::Dynamic, 1>
      hessian_diagonal = - n_samples_.cwiseProduct(elt_divide(exp_theta,
                                                   square(one + exp_theta)));
    hessian.resize(theta_size, theta_size);
    hessian.reserve(Eigen::VectorXi::Constant(theta_size, 1));
    for (int i = 0; i < theta_size; i++)
      hessian.insert(i, i) = hessian_diagonal(i);
  }

  /**
   * Returns the third derivative tensor. Because it is (cubic) diagonal,
   * the object is stored in a vector.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @param[in] eta_dummy additional likelihood parameters (used for other lk)
   * @return A vector containing the non-zero elements of the third
   *         derivative tensor.
   */
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta,
             const Eigen::Matrix<T2, Eigen::Dynamic, 1>& eta_dummy) const {
    Eigen::VectorXd exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());
    Eigen::VectorXd nominator = exp_theta.cwiseProduct(exp_theta - one);
    Eigen::VectorXd denominator = square(one + exp_theta)
                                   .cwiseProduct(one + exp_theta);

    return n_samples_.cwiseProduct(elt_divide(nominator, denominator));
  }

  Eigen::VectorXd compute_s2(const Eigen::VectorXd& theta,
                             const Eigen::VectorXd& eta,
                             const Eigen::MatrixXd& A,
                             int hessian_block_size) const {
    std::cout << "THIS FUNCTIONS SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::MatrixXd void_matrix;
    return void_matrix;
  }

  Eigen::VectorXd diff_eta_implicit(const Eigen::VectorXd& v,
                                    const Eigen::VectorXd& theta,
                                    const Eigen::VectorXd& eta) const {
    std::cout << "THIS FUNCTIONS SHOULD NEVER GET CALLED!" << std::endl;
    Eigen::MatrixXd void_matrix;
    return void_matrix;
  }
};

}  // namespace math
}  // namespace stan

#endif
