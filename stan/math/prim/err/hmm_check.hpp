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
 * @tparam T_omega type of the log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial guess vector
 * @param[in] log_omegas log matrix of observational densities.
 * @param[in] Gamma transition density between hidden states.
 * @param[in] rho initial state
 * @param[in] function the name of the function using the arguments.
 * @throw `std::invalid_argument` if Gamma is not square
 *         or if the size of rho is not the number of rows of log_omegas.
 * @throw `std::domain_error` if rho is not a simplex or the rows
 *         of Gamma are not a simplex.
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline void hmm_check(
    const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
    const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
    const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho, const char* function) {
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  check_consistent_size(function, "rho", rho, n_states);
  check_simplex(function, "rho", rho);
  check_square(function, "Gamma", Gamma);
  check_nonzero_size(function, "Gamma", Gamma);
  check_multiplicable(function, "Gamma", Gamma, "log_omegas", log_omegas);
  for (int i = 0; i < Gamma.rows(); ++i) {
    check_simplex(function, "Gamma[i, ]", row(Gamma, i + 1));
  }
}

}  // namespace math
}  // namespace stan
#endif
