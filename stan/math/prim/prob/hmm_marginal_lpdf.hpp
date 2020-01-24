#ifndef STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP
#define STAN_MATH_REV_FUN_HMM_MARGINAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <Eigen/Core>
#include <iostream>

namespace stan {
namespace math {
  // TO DO -- add checks for arguments.

  /**
   * For a Hidden Markov Model with observation y, hidden state x,
   * and parameters theta, return the log marginal density log
   * pi(y | theta). In this setting, the hidden states are discrete
   * and take values over the finite space {1, ..., K}.
   * The marginal lpdf is obtained via a forward pass.
   *
   * @param[in] log_omega log matrix of observational densities.
   *                      The (i, j)th entry corresponds to the
   *                      density of the ith observation, y_i,
   *                      should x_i = j.
   * @param[in] Gamma transition density between hidden states.
   *                  The (i, j)th entry is the probability x_n = i,
   *                  given x_{n - 1} = j.
   * @param[in] rho initial state
   * @param[in] n_states cardinality of the space for x.
   * @param[in, out] norm_norm log sum of the max coeff for each column
   *                           of alpha.
   * @return log marginal density.
   */
  double hmm_marginal_lpdf(const Eigen::MatrixXd& log_omegas,
                           const Eigen::MatrixXd& Gamma,
                           const Eigen::VectorXd& rho,
                           Eigen::MatrixXd& alphas,
                           Eigen::VectorXd& alpha_log_norms,
                           Eigen::MatrixXd& omegas) {
    omegas = log_omegas.array().exp();  // CHECK -- why the .array()?
    int n_transitions = log_omegas.rows() - 1;
    int n_states = log_omegas.cols();

    alphas.col(0) = omegas.col(0).cwiseProduct(rho);

    double norm = alphas.col(0).maxCoeff();
    alphas.col(0) /= norm;
    alpha_log_norms(0) = std::log(norm);

    for (int n = 0; n < n_transitions; n++) {
      alphas.col(n + 1)
        = omegas.col(n + 1).cwiseProduct(Gamma * alphas.col(n));

      double norm = alphas.col(n + 1).maxCoeff();
      alphas.col(n + 1) /= norm;
      alpha_log_norms(n + 1) = std::log(norm) + alpha_log_norms(n);
    }

    // CHECK -- do we need to "unnormalize" this?
    return alphas.col(n_transitions).sum();
  }


/**
 * Overload function for template types.
 */
template <typename T_omega, typename T_Gamma, typename T_rho>
inline return_type_t<T_omega, T_Gamma, T_rho> hmm_marginal_lpdf(
  const Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>& log_omegas,
  const Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>& Gamma,
  const Eigen::Matrix<T_rho, Eigen::Dynamic, 1>& rho) {
  // CHECK -- do I need to pass n_states ?
  using T_partials_return = partials_return_t<T_omega, T_Gamma, T_rho>;

  operands_and_partials<
    Eigen::Matrix<T_omega, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<T_Gamma, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Matrix<T_rho, Eigen::Dynamic, 1>
  > ops_partials(log_omegas, Gamma, rho);

  // operands_and_partials<T_omega, T_Gamma, T_rho>
  //   ops_partials(log_omegas, Gamma, rho);

  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;
  Eigen::MatrixXd alphas(n_states, n_transitions + 1);
  Eigen::VectorXd alpha_log_norms(n_transitions + 1);
  Eigen::MatrixXd omegas;

  T_partials_return log_marginal_density
    = hmm_marginal_lpdf(value_of(log_omegas),
                        value_of(Gamma),
                        value_of(rho),
                        alphas, alpha_log_norms, omegas);

  // Variables required for all three Jacobian-adjoint products.
  Eigen::VectorXd kappa = Eigen::VectorXd::Ones(n_states);
  Eigen::VectorXd kappa_log_norms(n_transitions);
  kappa_log_norms(n_transitions - 1) = 0;
  double norm_norm = alpha_log_norms(n_transitions);
  double grad_corr;

  if (!is_constant_all<T_Gamma>::value) {
    Eigen::MatrixXd Gamma_jacad(n_states, n_states);
    Gamma_jacad.setZero();
    grad_corr = std::exp(alpha_log_norms(n_transitions - 1) - norm_norm);

    Gamma_jacad += grad_corr
                     * kappa.cwiseProduct(omegas.col(n_transitions))
                     * alphas.col(n_transitions - 1).transpose();

    ops_partials.edge2_.partials_ = Gamma_jacad;
  }

  bool sensitivities_for_omega_or_rho
    = (!is_constant_all<T_omega>::value) || (!is_constant_all<T_rho>::value);

  // boundary terms
  if (sensitivities_for_omega_or_rho) {
    Eigen::MatrixXd log_omega_jacad(n_states, n_transitions + 1);
    log_omega_jacad.setZero();

    for (int n = n_transitions - 2; n >= 0; n--) {
      // CHECK -- is it worth creating a Gamma_dbl object?
      kappa = value_of(Gamma).transpose()
        * (omegas.col(n + 2).cwiseProduct(kappa));

      double norm = kappa.maxCoeff();
      kappa /= norm;
      kappa_log_norms(n) = std::log(norm) + kappa_log_norms(n + 1);

      if (!is_constant_all<T_omega>::value) {
        grad_corr = std::exp(alpha_log_norms(n) + kappa_log_norms(n + 1));

        log_omega_jacad.col(n + 1) = grad_corr
          * kappa.cwiseProduct(value_of(Gamma) * alphas.col(n));
      }
    }

    grad_corr = std::exp(kappa_log_norms(0) - norm_norm);
    Eigen::VectorXd c = value_of(Gamma).transpose()
      * omegas.col(1).cwiseProduct(kappa);

    if (!is_constant_all<T_omega>::value) {
      log_omega_jacad.col(0) = grad_corr * c.cwiseProduct(value_of(rho));
      ops_partials.edge1_.partials_ = log_omega_jacad;
    }

    if (!is_constant_all<T_rho>::value) {
      ops_partials.edge3_.partials_
        = grad_corr * c.cwiseProduct(omegas.col(0));
    }
  }

  return ops_partials.build(log_marginal_density);
}

}  // namespace math
}  // namespace stan
#endif
