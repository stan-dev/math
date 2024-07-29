#ifndef STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_BERNOULLI_LOGIT_LPMF_HPP
#define STAN_MATH_MIX_PROB_LAPLACE_MARGINAL_BERNOULLI_LOGIT_LPMF_HPP

#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>

namespace stan {
namespace math {

struct bernoulli_logit_likelihood {
  template <typename T_theta, typename T_eta>
  inline stan::return_type_t<T_theta, T_eta> operator()(
      const T_theta& theta, const T_eta& /* eta */, const Eigen::VectorXd& y,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    return sum(elt_multiply(theta, y)
               - elt_multiply(to_vector(delta_int), log(add(1.0, exp(theta)))));
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
  laplace_options ops{hessian_block_size, solver, max_steps_line_search,
                      tolerance, max_num_steps};
  return laplace_marginal_density(bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples),
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
  constexpr laplace_options ops{1, 1, 0, 1e-6, 100};
  return laplace_marginal_density(bernoulli_logit_likelihood{},
      std::forward_as_tuple(to_vector(y), n_samples),
      covariance_function, eta_dummy, theta_0, msgs, ops,
      std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
