#ifndef STAN_MATH_PRIM_PROB_HMM_LATENT_RNG_HPP
#define STAN_MATH_PRIM_PROB_HMM_LATENT_RNG_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/hmm_check.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/random.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * For a hidden Markov model with observation y, hidden state x,
 * and parameters theta, generate samples from the posterior distribution
 * of the hidden states, x.
 * In this setting, the hidden states are discrete
 * and takes values over the finite space {1, ..., K}.
 * log_omegas is a matrix of observational densities, where
 * the (i, j)th entry corresponds to the density of the ith observation, y_i,
 * given x_i = j.
 * The transition matrix Gamma is such that the (i, j)th entry is the
 * probability that x_n = j given x_{n - 1} = i. The rows of Gamma are
 * simplexes.
 *
 * @tparam T_omega type of the log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial guess vector
 * @param[in] log_omegas log matrix of observational densities.
 * @param[in] Gamma transition density between hidden states.
 * @param[in] rho initial state
 * @param[in] rng random number generator
 * @return sample from the posterior distribution of the hidden states.
 * @throw `std::invalid_argument` if Gamma is not square, when we have
 *         at least one transition, or if the size of rho is not the
 *         number of rows of log_omegas.
 * @throw `std::domain_error` if rho is not a simplex and of the rows
 *         of Gamma are not a simplex (when there is at least one transition).
 */
template <typename T_omega, typename T_Gamma, typename T_rho, class RNG,
          require_all_eigen_t<T_omega, T_Gamma>* = nullptr,
          require_eigen_col_vector_t<T_rho>* = nullptr>
inline std::vector<int> hmm_latent_rng(const T_omega& log_omegas,
                                       const T_Gamma& Gamma, const T_rho& rho,
                                       RNG& rng) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  Eigen::MatrixXd omegas = value_of(log_omegas).array().exp();
  ref_type_t<decltype(value_of(rho))> rho_dbl = value_of(rho);
  ref_type_t<decltype(value_of(Gamma))> Gamma_dbl = value_of(Gamma);
  hmm_check(log_omegas, Gamma_dbl, rho_dbl, "hmm_latent_rng");

  Eigen::MatrixXd alphas(n_states, n_transitions + 1);
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_dbl);
  alphas.col(0) /= alphas.col(0).maxCoeff();

  for (int n = 0; n < n_transitions; ++n) {
    alphas.col(n + 1)
        = omegas.col(n + 1).cwiseProduct(Gamma_dbl.transpose() * alphas.col(n));
    alphas.col(n + 1) /= alphas.col(n + 1).maxCoeff();
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Ones(n_states);

  // sample the last hidden state
  std::vector<int> hidden_states(n_transitions + 1);
  std::vector<double> probs(n_states);
  Eigen::Map<Eigen::VectorXd> probs_vec(probs.data(), n_states);
  probs_vec = alphas.col(n_transitions) / alphas.col(n_transitions).sum();
  boost::random::discrete_distribution<> cat_hidden(probs);
  hidden_states[n_transitions] = cat_hidden(rng) + stan::error_index::value;

  for (int n = n_transitions; n-- > 0;) {
    // Sample the nth hidden state conditional on (n + 1)st hidden state.
    // Subtract error_index in order to use C++ index.
    int last_hs = hidden_states[n + 1] - stan::error_index::value;

    probs_vec = alphas.col(n).cwiseProduct(Gamma_dbl.col(last_hs))
                * beta(last_hs) * omegas(last_hs, n + 1);

    probs_vec /= probs_vec.sum();

    // discrete_distribution produces samples in [0, K), so
    // we need to add 1 to generate over [1, K).
    boost::random::discrete_distribution<> cat_hidden(probs);
    hidden_states[n] = cat_hidden(rng) + stan::error_index::value;

    // update backwards state
    beta = Gamma_dbl * (omegas.col(n + 1).cwiseProduct(beta));
    beta /= beta.maxCoeff();
  }

  return hidden_states;
}

}  // namespace math
}  // namespace stan
#endif
