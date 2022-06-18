#ifndef STAN_MATH_PRIM_ERR_HMM_CHECK_HPP
#define STAN_MATH_PRIM_ERR_HMM_CHECK_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

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
template <typename T_omega, typename T_Gamma, typename T_rho,
          require_all_eigen_t<T_omega, T_Gamma>* = nullptr,
          require_eigen_col_vector_t<T_rho>* = nullptr>
inline void hmm_check(const T_omega& log_omegas, const T_Gamma& Gamma,
                      const T_rho& rho, const char* function) {
  int n_states = log_omegas.rows();

  check_consistent_size(function, "rho", rho, n_states);
  check_square(function, "Gamma", Gamma);
  check_nonzero_size(function, "Gamma", Gamma);
  check_multiplicable(function, "Gamma", Gamma, "log_omegas", log_omegas);

  check_simplex(function, "rho", rho);
  const auto& Gamma_ref = to_ref(Gamma);
  for (int i = 0; i < Gamma.rows(); ++i) {
    check_simplex(function, "Gamma[i, ]", Gamma_ref.row(i));
  }
}

}  // namespace math
}  // namespace stan
#endif
