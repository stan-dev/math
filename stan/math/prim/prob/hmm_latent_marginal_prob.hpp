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
 * and parameters theta, compute the marginal probability
 * of the hidden state.
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
  alphas.col(0) = omegas.col(0).cwiseProduct(rho);
  alphas.col(0) /= alphas.col(0).maxCoeff();

  for (int n = 0; n < n_transitions; ++n)
    alphas.col(n + 1) = omegas.col(n + 1).cwiseProduct(Gamma * alphas.col(n));

  // Backward pass with running normalization
  Eigen::VectorXd beta = Eigen::VectorXd::Ones(n_states);

  alphas.col(n_transitions) /= alphas.col(n_transitions).sum();

  for (int n = n_transitions - 1; n >= 0; --n) {
    beta = Gamma.transpose() * (omegas.col(n + 1).cwiseProduct(beta));
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
