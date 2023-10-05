#ifndef STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP
#define STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/row.hpp>
#include <stan/math/prim/fun/col.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T_omega, typename T_Gamma, typename T_rho, typename T_alpha>
inline auto hmm_marginal_val(
    const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& omegas,
    const T_Gamma& Gamma_val, const T_rho& rho_val,
    Eigen::Matrix<T_alpha, Eigen::Dynamic, Eigen::Dynamic>& alphas,
    Eigen::Matrix<T_alpha, Eigen::Dynamic, 1>& alpha_log_norms,
    T_alpha& norm_norm) {
  const int n_transitions = omegas.cols() - 1;
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_val);

  const auto norm = alphas.col(0).maxCoeff();
  alphas.col(0) /= norm;
  alpha_log_norms(0) = log(norm);

  auto Gamma_val_transpose = Gamma_val.transpose().eval();
  for (int n = 1; n <= n_transitions; ++n) {
    alphas.col(n)
        = omegas.col(n).cwiseProduct(Gamma_val_transpose * alphas.col(n - 1));
    const auto col_norm = alphas.col(n).maxCoeff();
    alphas.col(n) /= col_norm;
    alpha_log_norms(n) = log(col_norm) + alpha_log_norms(n - 1);
  }
  norm_norm = alpha_log_norms(n_transitions);
  return log(alphas.col(n_transitions).sum()) + norm_norm;
}

/**
 * For a Hidden Markov Model with observation y, hidden state x,
 * and parameters theta, return the log marginal density, log
 * p(y | theta). In this setting, the hidden states are discrete
 * and take values over the finite space {1, ..., K}.
 * The marginal lpdf is obtained via a forward pass, and
 * the derivative is calculated with an adjoint method,
 * e.g (Betancourt, Margossian, & Leos-Barajas, 2020).
 * log_omegas is a matrix of observational densities, where
 * the (i, j)th entry corresponds to the density of the jth observation, y_j,
 * given x_j = i.
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
 * @return log marginal density.
 * @throw `std::invalid_argument` if Gamma is not square, when we have
 *         at least one transition, or if the size of rho is not the
 *         number of rows of log_omegas.
 * @throw `std::domain_error` if rho is not a simplex and of the rows
 *         of Gamma are not a simplex (when there is at least one transition).
 */
template <typename T_omega, typename T_Gamma, typename T_rho,
          require_all_eigen_t<T_omega, T_Gamma>* = nullptr,
          require_eigen_col_vector_t<T_rho>* = nullptr>
inline auto hmm_marginal(const T_omega& log_omegas, const T_Gamma& Gamma,
                         const T_rho& rho) {
  using T_partial_type = partials_return_t<T_omega, T_Gamma, T_rho>;
  using eig_matrix_partial
      = Eigen::Matrix<T_partial_type, Eigen::Dynamic, Eigen::Dynamic>;
  using eig_vector_partial = Eigen::Matrix<T_partial_type, Eigen::Dynamic, 1>;
  using T_omega_ref = ref_type_if_not_constant_t<T_omega>;
  using T_Gamma_ref = ref_type_if_not_constant_t<T_Gamma>;
  using T_rho_ref = ref_type_if_not_constant_t<T_rho>;
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  T_omega_ref log_omegas_ref = log_omegas;
  T_Gamma_ref Gamma_ref = Gamma;
  T_rho_ref rho_ref = rho;

  const auto& Gamma_val = to_ref(value_of(Gamma_ref));
  const auto& rho_val = to_ref(value_of(rho_ref));
  hmm_check(log_omegas, Gamma_val, rho_val, "hmm_marginal");

  auto ops_partials
      = make_partials_propagator(log_omegas_ref, Gamma_ref, rho_ref);

  eig_matrix_partial alphas(n_states, n_transitions + 1);
  eig_vector_partial alpha_log_norms(n_transitions + 1);
  // compute the density using the forward algorithm.
  eig_matrix_partial omegas = value_of(log_omegas_ref).array().exp();
  T_partial_type norm_norm;
  auto log_marginal_density = hmm_marginal_val(
      omegas, Gamma_val, rho_val, alphas, alpha_log_norms, norm_norm);

  // Variables required for all three Jacobian-adjoint products.
  auto unnormed_marginal = alphas.col(n_transitions).sum();

  std::vector<eig_vector_partial> kappa(n_transitions);
  eig_vector_partial kappa_log_norms(n_transitions);
  std::vector<T_partial_type> grad_corr(n_transitions, 0);

  if (n_transitions > 0) {
    kappa[n_transitions - 1] = Eigen::VectorXd::Ones(n_states);
    kappa_log_norms(n_transitions - 1) = 0;
    grad_corr[n_transitions - 1]
        = exp(alpha_log_norms(n_transitions - 1) - norm_norm);
  }

  for (int n = n_transitions - 1; n-- > 0;) {
    kappa[n] = Gamma_val * (omegas.col(n + 2).cwiseProduct(kappa[n + 1]));

    auto norm = kappa[n].maxCoeff();
    kappa[n] /= norm;
    kappa_log_norms[n] = log(norm) + kappa_log_norms[n + 1];
    grad_corr[n] = exp(alpha_log_norms[n] + kappa_log_norms[n] - norm_norm);
  }

  if (!is_constant_all<T_Gamma>::value) {
    for (int n = n_transitions - 1; n >= 0; --n) {
      edge<1>(ops_partials).partials_
          += grad_corr[n] * alphas.col(n)
             * kappa[n].cwiseProduct(omegas.col(n + 1)).transpose()
             / unnormed_marginal;
    }
  }

  if (!is_constant_all<T_omega, T_rho>::value) {
    // Boundary terms
    if (n_transitions == 0) {
      if (!is_constant_all<T_omega>::value) {
        edge<0>(ops_partials).partials_
            = omegas.cwiseProduct(rho_val) / exp(log_marginal_density);
      }

      if (!is_constant_all<T_rho>::value) {
        edge<2>(ops_partials).partials_
            = omegas.col(0) / exp(log_marginal_density);
      }
      return ops_partials.build(log_marginal_density);
    } else {
      auto grad_corr_boundary = exp(kappa_log_norms(0) - norm_norm);
      eig_vector_partial C = Gamma_val * omegas.col(1).cwiseProduct(kappa[0]);

      if (!is_constant_all<T_omega>::value) {
        eig_matrix_partial log_omega_jacad
            = Eigen::MatrixXd::Zero(n_states, n_transitions + 1);

        for (int n = n_transitions - 1; n >= 0; --n) {
          log_omega_jacad.col(n + 1)
              = grad_corr[n]
                * kappa[n].cwiseProduct(Gamma_val.transpose() * alphas.col(n));
        }

        log_omega_jacad.col(0) = grad_corr_boundary * C.cwiseProduct(rho_val);
        edge<0>(ops_partials).partials_
            = log_omega_jacad.cwiseProduct(omegas / unnormed_marginal);
      }

      if (!is_constant_all<T_rho>::value) {
        partials<2>(ops_partials) = grad_corr_boundary
                                    * C.cwiseProduct(omegas.col(0))
                                    / unnormed_marginal;
      }
    }
  }

  return ops_partials.build(log_marginal_density);
}

}  // namespace math
}  // namespace stan
#endif
