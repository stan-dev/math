#ifndef STAN_MATH_PRIM_PROB_HMM_LATENT_MARGINAL_PROB_HPP
#define STAN_MATH_PRIM_PROB_HMM_LATENT_MARGINAL_PROB_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/hmm_check.hpp>
#include <Eigen/Core>
#include <boost/random.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * For a hidden Markov model with observation y, hidden state x,
 * and parameters theta, compute the marginal posterior probability
 * of the hidden states. In this setting, the hidden states are discrete
 * and takes values over the finite space {1, ..., K}.
 * log_omegas is a matrix of observational densities, where
 * the (i, j)th entry corresponds to the density of the ith observation, y_i,
 * given x_i = j.
 * The transition matrix Gamma is such that the (i, j)th entry is the
 * probability that x_n = j given x_{n - 1} = i. The rows of Gamma are
 * simplexes.
 * Note we only calculate marginal probabilities. However correlation between
 * the latent states must be considered in order to sample latent states.
 *
 * @tparam T_omega type of the log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial guess vector
 * @param[in] log_omegas log matrix of observational densities.
 * @param[in] Gamma transition density between hidden states.
 * @param[in] rho initial state
 * @return the posterior probability for each latent state.
 * @throw `std::invalid_argument` if Gamma is not square, when we have
 *         at least one transition, or if the size of rho is not the
 *         number of rows of log_omegas.
 * @throw `std::domain_error` if rho is not a simplex and of the rows
 *         of Gamma are not a simplex (when there is at least one transition).
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline Eigen::MatrixXd hmm_latent_marginal_prob(
  const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
  const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
  const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  hmm_check(log_omegas, Gamma, rho, "hmm_latent_marginal_prob");

  Eigen::MatrixXd omegas = value_of(log_omegas).array().exp();
  Eigen::VectorXd rho_dbl = value_of(rho);
  Eigen::MatrixXd Gamma_dbl = value_of(Gamma);

  Eigen::MatrixXd alphas(n_states, n_transitions + 1);
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_dbl);
  alphas.col(0) /= alphas.col(0).maxCoeff();

  Eigen::MatrixXd Gamma_dbl_transpose = Gamma_dbl.transpose();
  for (int n = 1; n <= n_transitions; ++n)
    alphas.col(n) = omegas.col(n).cwiseProduct(Gamma_dbl_transpose
                      * alphas.col(n - 1));

  // Backward pass with running normalization
  Eigen::VectorXd beta = Eigen::VectorXd::Ones(n_states);

  alphas.col(n_transitions) /= alphas.col(n_transitions).sum();

  for (int n = n_transitions - 1; n >= 0; --n) {
    beta = Gamma_dbl * (omegas.col(n + 1).cwiseProduct(beta));
    beta /= beta.maxCoeff();

    // Reuse alphas to store probabilities
    alphas.col(n) = alphas.col(n).cwiseProduct(beta);
    alphas.col(n) /= alphas.col(n).sum();
  }

  return alphas;
}

}  // namespace math
}  // namespace stan
#endif
