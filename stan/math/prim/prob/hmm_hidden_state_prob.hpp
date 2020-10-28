#ifndef STAN_MATH_PRIM_PROB_HMM_HIDDEN_STATE_PROB_HPP
#define STAN_MATH_PRIM_PROB_HMM_HIDDEN_STATE_PROB_HPP

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
 * and parameters theta, compute the marginal probability
 * vector for each x, given y and theta, p(x_i | y, theta).
 * In this setting, the hidden states are discrete
 * and take values over the finite space {1, ..., K}.
 * Hence for each hidden variable x, we compute a simplex with K elements.
 * The final result is stored in a K by N matrix, where N is the length of x.
 * log_omegas is a matrix of observational densities, where
 * the (i, j)th entry corresponds to the density of the ith observation, y_i,
 * given x_i = j.
 * The transition matrix Gamma is such that the (i, j)th entry is the
 * probability that x_n = j given x_{n - 1} = i. The rows of Gamma are
 * simplexes.
 * This function cannot be used to reconstruct the marginal distributon
 * of a state sequence given parameters and an observation sequence,
 * p(x | y, theta),
 * because it only computes marginals on a state-by-state basis.
 *
 * @tparam T_omega type of the log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial guess vector
 * @param[in] log_omegas log matrix of observational densities
 * @param[in] Gamma transition density between hidden states
 * @param[in] rho initial state
 * @return the posterior probability for each latent state
 * @throw `std::invalid_argument` if Gamma is not square
 *         or if the size of rho is not the number of rows of log_omegas
 * @throw `std::domain_error` if rho is not a simplex and of the rows
 *         of Gamma are not a simplex
 */
template <typename T_omega, typename T_Gamma, typename T_rho,
          require_all_eigen_t<T_omega, T_Gamma>* = nullptr,
          require_eigen_col_vector_t<T_rho>* = nullptr>
inline Eigen::MatrixXd hmm_hidden_state_prob(const T_omega& log_omegas,
                                             const T_Gamma& Gamma,
                                             const T_rho& rho) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  Eigen::MatrixXd omegas = value_of(log_omegas).array().exp();
  ref_type_t<decltype(value_of(rho))> rho_dbl = value_of(rho);
  ref_type_t<decltype(value_of(Gamma))> Gamma_dbl = value_of(Gamma);
  hmm_check(log_omegas, Gamma_dbl, rho_dbl, "hmm_hidden_state_prob");

  Eigen::MatrixXd alphas(n_states, n_transitions + 1);
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_dbl);
  alphas.col(0) /= alphas.col(0).maxCoeff();

  for (int n = 0; n < n_transitions; ++n)
    alphas.col(n + 1)
        = omegas.col(n + 1).cwiseProduct(Gamma_dbl.transpose() * alphas.col(n));

  // Backward pass with running normalization
  Eigen::VectorXd beta = Eigen::VectorXd::Ones(n_states);

  alphas.col(n_transitions) /= alphas.col(n_transitions).sum();

  for (int n = n_transitions; n-- > 0;) {
    beta = Gamma_dbl * omegas.col(n + 1).cwiseProduct(beta);
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
