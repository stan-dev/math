#ifndef STAN_MATH_PRIM_ERR_HMM_CHECK_HPP
#define STAN_MATH_PRIM_ERR_HMM_CHECK_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/row.hpp>

namespace stan {
namespace math {

/**
 * Check arguments for hidden Markov model functions with a discrete
 * latent state (lpdf, rng for latent states, and marginal probabilities
 * for latent sates).
 *
 * @tparam T_omega type of log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial state vector.
 *
 * @param[in] log_omega log matrix of observational densities.
 *              The (i, j)th entry corresponds to the
 *              density of the ith observation, y_i,
 *              given x_i = j.
 * @param[in] Gamma transition density between hidden states.
 *              The (i, j)th entry is the probability that x_n = j,
 *              given x_{n - 1} = i. The rows of Gamma are simplexes.
 * @param[in] rho initial state
 * @throw `std::domain_error` if Gamma is not square.
 * @throw `std::invalid_argument` if the size of rho is not
 * the number of rows of Gamma.
 * @throw `std::domain_error` if rho is not a simplex.
 * @throw `std::domain_error` if the rows of Gamma are
 * not a simplex.
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline void hmm_check(
  const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
  const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
  const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  check_square("hmm_marginal_lpdf", "Gamma", Gamma);
  check_consistent_size("hmm_marginal_lpdf", "Gamma", row(Gamma, 1), n_states);
  check_consistent_size("hmm_marginal_lpdf", "rho", rho, n_states);
  check_simplex("hmm_marginal_lpdf", "rho", rho);
  for (int i = 0; i < Gamma.rows(); ++i) {
    check_simplex("hmm_marginal_lpdf", "Gamma[i, ]", row(Gamma, i + 1));
  }
}

}  // namespace math
}  // namespace stan
#endif
