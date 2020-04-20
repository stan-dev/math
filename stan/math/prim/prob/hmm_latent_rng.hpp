#ifndef STAN_MATH_PRIM_PROB_HMM_LATENT_RNG_HPP
#define STAN_MATH_PRIM_PROB_HMM_LATENT_RNG_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <Eigen/Core>
#include <boost/random.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * For a hidden Markov model with observation y, hidden state x,
 * and parameters theta, sample hidden states x.
 */
template <typename T_omega, typename T_Gamma, typename T_rho, class RNG>
inline std::vector<int> hmm_latent_rng(
  const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
  const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
  const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho,
  RNG& rng) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  // TODO (charlesm93): Add checks. Create an hmm_check function
  // to be used for all hmm function.

  Eigen::MatrixXd omegas = value_of(log_omegas).array().exp();
  Eigen::VectorXd rho_dbl = value_of(rho);
  Eigen::MatrixXd Gamma_dbl = value_of(Gamma);

  Eigen::MatrixXd alphas(n_states, n_transitions + 1);
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_dbl);
  alphas.col(0) /= alphas.col(0).maxCoeff();

  for (int n = 0; n < n_transitions; ++n) {
    alphas.col(n + 1) = omegas.col(n + 1).
                          cwiseProduct(Gamma_dbl * alphas.col(n));
    alphas.col(n + 1) /= alphas.col(n + 1).maxCoeff();
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Ones(n_states);

  // sample the last hidden state
  std::vector<int> hidden_states(n_transitions + 1);
  Eigen::VectorXd probs_vec = alphas.col(n_transitions)
                              / alphas.col(n_transitions).sum();
  std::vector<double> probs(probs_vec.data(), probs_vec.data() + n_states);
  boost::random::discrete_distribution<> cat_hidden(probs);
  hidden_states[n_transitions] = cat_hidden(rng);

  for (int n = n_transitions - 1; n >= 0; --n) {
    // sample the nth hidden state conditional on (n + 1)st hidden state
    int last_hs = hidden_states[n + 1];

    // CHECK -- do we need to redeclare the hidden state?
    for (int k = 0; k < n_states; ++k) {
      probs_vec[k] = alphas(k, n) * omegas(last_hs, n + 1)
                     * Gamma_dbl(last_hs, k) * beta(last_hs);
    }
    probs_vec /= probs_vec.sum();
    std::vector<double> probs(probs_vec.data(), probs_vec.data() + n_states);
    boost::random::discrete_distribution<> cat_hidden(probs);
    hidden_states[n] = cat_hidden(rng);

    // update backwards state
    beta = Gamma.transpose() * (omegas.col(n + 1).cwiseProduct(beta));
    beta /= beta.maxCoeff();
  }

  return hidden_states;
}

}  // namespace math
}  // namespace stan
#endif
