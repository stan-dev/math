#ifndef STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP
#define STAN_MATH_LAPLACE_LAPLACE_MARGINAL_BERNOULLI_LPMF_HPP

#include <stan/math/mix/laplace/laplace_marginal.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct bernoulli_logit_likelihood {
  template <typename T_theta, typename T_eta>
  inline stan::return_type_t<T_theta, T_eta> operator()(
      const T_theta& theta, const T_eta& /* eta */, const Eigen::VectorXd& y,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    return sum(theta.cwiseProduct(y)
               - to_vector(delta_int).cwiseProduct(log(add(1.0, exp(theta)))));
  }
};

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
  template <typename Vec1, typename Vec2>
  diff_bernoulli_logit(Vec1&& n_samples, Vec2&& sums)
      : n_samples_(std::forward<Vec1>(n_samples)),
        sums_(std::forward<Vec2>(sums)) {}

  /**
   * Return the log density.
   * @tparam T type of the log poisson parameter.
   * @param[in] theta log poisson parameters for each group.
   * @return the log density.
   */
  template <typename Theta, typename Eta,
            require_eigen_vector_t<Theta>* = nullptr,
            require_eigen_t<Eta>* = nullptr>
  inline auto log_likelihood(const Theta& theta,
                             const Eta& /* eta_dummy */) const {
    return sum(theta.cwiseProduct(sums_)
               - n_samples_.cwiseProduct(
                   log(add(Eigen::VectorXd::Ones(theta.size()), exp(theta)))));
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
  template <typename Theta, typename Eta, typename GradVec,
            require_eigen_vector_t<Theta>* = nullptr,
            require_eigen_t<Eta>* = nullptr>
  inline Eigen::SparseMatrix<double> diff(const Theta& theta,
                                          const Eta& /* eta */,
                                          GradVec& gradient,
                                          const Eigen::Index block_size_dummy
                                          = 0) const {
    auto exp_theta = exp(theta).eval();
    const Eigen::Index theta_size = theta.size();
    gradient = sums_ - n_samples_.cwiseProduct(inv_logit(theta));
    auto hessian_diagonal = (-n_samples_.cwiseProduct(elt_divide(
                                 exp_theta, square(add(1.0, exp_theta)))))
                                .eval();
    Eigen::SparseMatrix<double> hessian(theta_size, theta_size);
    hessian.reserve(Eigen::VectorXi::Constant(theta_size, 1));
    for (Eigen::Index i = 0; i < theta_size; i++) {
      hessian.insert(i, i) = hessian_diagonal(i);
    }
    return hessian;
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
  template <typename Theta, typename Eta,
            require_eigen_vector_t<Theta>* = nullptr,
            require_eigen_t<Eta>* = nullptr>
  inline auto third_diff(const Theta& theta, const Eta& /* eta */) const {
    Eigen::VectorXd exp_theta = exp(theta);
    auto nominator = exp_theta.cwiseProduct(
        exp_theta - Eigen::VectorXd::Ones(theta.size()));
    Eigen::VectorXd one_plus_exp_theta = add(1.0, exp_theta);
    auto denominator
        = square(one_plus_exp_theta).cwiseProduct(one_plus_exp_theta);

    return (n_samples_.cwiseProduct(elt_divide(nominator, denominator))).eval();
  }
};

/**
 * Wrapper function around the laplace_marginal function for
 * a logistic Bernoulli likelihood. Returns the marginal density
 * p(y | phi) by marginalizing out the latent gaussian variable,
 * with a Laplace approximation. See the laplace_marginal function
 * for more details.
 *
 * @tparam T0 The type of the initial guess, theta_0.
 * @tparam T1 The type for the global parameter, phi.
 * @param[in] theta_0 the initial guess for the Laplace approximation.
 * @param[in] phi model parameters for the covariance function.
 * @param[in] x data for the covariance function.
 * @param[in] n_samples number of samples per group. First sufficient
 *            statistics.
 * @param[in] y total counts per group. Second sufficient statistics.
 * @param[in] tolerance controls the convergence criterion when finding
 *            the mode in the Laplace approximation.
 * @param[in] max_num_steps maximum number of steps before the Newton solver
 *            breaks and returns an error.
 */
template <typename CovarF, typename ThetaMatrix, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline auto laplace_marginal_tol_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    double tolerance, long int max_num_steps, const int hessian_block_size,
    const int solver, const int max_steps_line_search,
    const ThetaMatrix& theta_0, CovarF&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  laplace_options ops{hessian_block_size, solver,
    max_steps_line_search, tolerance, max_num_steps};
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

template <typename CovarF, typename ThetaMatrix, typename... Args,
          require_eigen_t<ThetaMatrix>* = nullptr>
inline auto laplace_marginal_bernoulli_logit_lpmf(
    const std::vector<int>& y, const std::vector<int>& n_samples,
    const ThetaMatrix& theta_0, CovarF&& covariance_function,
    std::ostream* msgs, Args&&... args) {
  // TODO: change this to a VectorXd once we have operands & partials.
  Eigen::Matrix<double, 0, 0> eta_dummy;
  laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(
      diff_likelihood<bernoulli_logit_likelihood>(
          bernoulli_logit_likelihood{}, to_vector(y), n_samples, msgs),
      covariance_function, eta_dummy, theta_0, msgs,
      std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
