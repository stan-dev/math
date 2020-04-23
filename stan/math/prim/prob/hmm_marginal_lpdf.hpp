#ifndef STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP
#define STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/row.hpp>
#include <stan/math/prim/fun/col.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/core.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T_omega, typename T_Gamma, typename T_rho, typename T_alpha>
inline auto hmm_marginal_lpdf_val(
    const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& omegas,
    const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma_val,
    const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho_val,
    Eigen::Matrix<T_alpha, Eigen::Dynamic, Eigen::Dynamic>& alphas,
    Eigen::Matrix<T_alpha, Eigen::Dynamic, 1>& alpha_log_norms) {
  const int n_states = omegas.rows();
  const int n_transitions = omegas.cols() - 1;
  alphas.col(0) = omegas.col(0).cwiseProduct(rho_val);

  const auto norm = alphas.col(0).maxCoeff();
  alphas.col(0) /= norm;
  alpha_log_norms(0) = log(norm);

  for (int n = 1; n <= n_transitions; ++n) {
    alphas.col(n)
        = omegas.col(n).cwiseProduct(Gamma_val.transpose() * alphas.col(n - 1));
    const auto col_norm = alphas.col(n).maxCoeff();
    alphas.col(n) /= col_norm;
    alpha_log_norms(n) = log(col_norm) + alpha_log_norms(n - 1);
  }
  return log(alphas.col(n_transitions).sum()) + alpha_log_norms(n_transitions);
}

/**
 * For a Hidden Markov Model with observation y, hidden state x,
 * and parameters theta, return the log marginal density, log
 * p(y | theta). In this setting, the hidden states are discrete
 * and take values over the finite space {1, ..., K}.
 * The marginal lpdf is obtained via a forward pass, and
 * the derivative is calculated with an adjoint method,
 * see (Betancourt, Margossian, & Leos-Barajas, 2020).
 * The Gamma argument is only checked if there is at least one transition.
 *
 * @tparam T_omega type of the log likelihood matrix
 * @tparam T_Gamma type of the transition matrix
 * @tparam T_rho type of the initial guess vector
 * @param[in] log_omegas log matrix of observational densities.
 *              The (i, j)th entry corresponds to the
 *              density of the ith observation, y_i,
 *              given x_i = j.
 * @param[in] Gamma transition density between hidden states.
 *              The (i, j)th entry is the probability that x_n = j,
 *              given x_{n - 1} = i. The rows of Gamma are simplexes.
 * @param[in] rho initial state
 * @return log marginal density.
 * @throw `std::invalid_argument` if Gamma is not square, when we have
 *         at least one transition.
 * @throw `std::invalid_argument` if the size of rho is not
 * the number of rows of log_omegas.
 * @throw `std::domain_error` if rho is not a simplex and of the rows
 *         of Gamma are not a simplex (when there is at least one transition).
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline auto hmm_marginal_lpdf(
    const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
    const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
    const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho) {
  using T_partial_type = partials_return_t<T_omega, T_Gamma, T_rho>;
  using eig_matrix_partial
      = Eigen::Matrix<T_partial_type, Eigen::Dynamic, Eigen::Dynamic>;
  using eig_vector_partial = Eigen::Matrix<T_partial_type, Eigen::Dynamic, 1>;
  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  check_consistent_size("hmm_marginal_lpdf", "rho", rho, n_states);
  check_simplex("hmm_marginal_lpdf", "rho", rho);
  if (n_transitions != 0) {
    check_square("hmm_marginal_lpdf", "Gamma", Gamma);
    check_nonzero_size("hmm_marginal_lpdf", "Gamma", Gamma);
    check_matching_sizes("hmm_marginal_lpdf", "Gamma (row and column)",
                         Gamma.row(0), "log_omegas (row)", log_omegas.col(0));
    for (int i = 0; i < Gamma.rows(); ++i) {
      check_simplex("hmm_marginal_lpdf", "Gamma[i, ]", row(Gamma, i + 1));
    }
  }

  operands_and_partials<Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>,
                        Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>,
                        Eigen::Matrix<T_rho, Eigen::Dynamic, 1> >
      ops_partials(log_omegas, Gamma, rho);

  eig_matrix_partial alphas(n_states, n_transitions + 1);
  eig_vector_partial alpha_log_norms(n_transitions + 1);
  auto Gamma_val = value_of(Gamma);

  // compute the density using the forward algorithm.
  auto rho_val = value_of(rho);
  eig_matrix_partial omegas = value_of(log_omegas).array().exp();
  auto log_marginal_density = hmm_marginal_lpdf_val(omegas, Gamma_val, rho_val,
                                                    alphas, alpha_log_norms);

  // Variables required for all three Jacobian-adjoint products.
  auto norm_norm = alpha_log_norms(n_transitions);
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
    eig_matrix_partial Gamma_jacad = Eigen::MatrixXd::Zero(n_states, n_states);

    for (int n = n_transitions - 1; n >= 0; --n) {
      Gamma_jacad += grad_corr[n] * alphas.col(n)
                     * kappa[n].cwiseProduct(omegas.col(n + 1)).transpose()
                     / unnormed_marginal;
    }

    ops_partials.edge2_.partials_ = Gamma_jacad;
  }

  if (!is_constant_all<T_omega, T_rho>::value) {
    eig_matrix_partial log_omega_jacad
        = Eigen::MatrixXd::Zero(n_states, n_transitions + 1);

    if (!is_constant_all<T_omega>::value) {
      // auto Gamma_alpha = Gamma_val.transpose() * alphas;
      for (int n = n_transitions - 1; n >= 0; --n)
        log_omega_jacad.col(n + 1)
            = grad_corr[n]
              // * kappa[n].cwiseProduct(Gamma_alpha.col(n)).eval();
              //  * kappa[n].cwiseProduct(Gamma_alpha.col(n));
              * kappa[n].cwiseProduct(Gamma_val.transpose() * alphas.col(n));
    }

    // Boundary terms
    if (n_transitions == 0) {
      if (!is_constant_all<T_omega>::value) {
        log_omega_jacad
            = omegas.cwiseProduct(value_of(rho)) / exp(log_marginal_density);
        ops_partials.edge1_.partials_ = log_omega_jacad;
      }

      if (!is_constant_all<T_rho>::value) {
        ops_partials.edge3_.partials_
            = omegas.col(0) / exp(log_marginal_density);
      }
      return ops_partials.build(log_marginal_density);
    } else {
      auto grad_corr_boundary = exp(kappa_log_norms(0) - norm_norm);
      eig_vector_partial C = Gamma_val * omegas.col(1).cwiseProduct(kappa[0]);

      if (!is_constant_all<T_omega>::value) {
        log_omega_jacad.col(0)
            = grad_corr_boundary * C.cwiseProduct(value_of(rho));
        ops_partials.edge1_.partials_
            = log_omega_jacad.cwiseProduct(omegas / unnormed_marginal);
      }

      if (!is_constant_all<T_rho>::value) {
        ops_partials.edge3_.partials_ = grad_corr_boundary
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
