#ifndef STAN_MATH_OPENCL_PRIM_HMM_MARGINAL_HPP
#define STAN_MATH_OPENCL_PRIM_HMM_MARGINAL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
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
template <
    bool propto, typename T_omega_cl, typename T_Gamma_cl, typename T_rho_cl,
    require_all_prim_or_rev_kernel_expression_t<T_omega_cl, T_Gamma_cl,
                                                T_rho_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_omega_cl, T_Gamma_cl, T_rho_cl>* = nullptr>
inline return_type_t<T_omega_cl, T_Gamma_cl, T_rho_cl> hmm_marginal(
    const T_omega_cl& log_omegas, const T_Gamma_cl& Gamma,
    const T_rho_cl& rho) {
  static const char* function = "hmm_marginal(OpenCL)";
  using T_partials_return = partials_return_t<T_omega_cl, T_Gamma_cl, T_rho_cl>;
  using std::isfinite;
  using std::isnan;

  int n_states = log_omegas.rows();
  int n_transitions = log_omegas.cols() - 1;

  check_vector(function, "rho", rho);
  check_consistent_size(function, "rho", rho, n_states);
  check_square(function, "Gamma", Gamma);
  check_nonzero_size(function, "Gamma", Gamma);
  check_multiplicable(function, "Gamma", Gamma, "log_omegas", log_omegas);

  // rho simplex
  // gamma rows simlex

  const auto& log_omegas_val = value_of(log_omegas);
  const auto& Gamma_val = value_of(Gamma);
  const auto& rho_val = value_of(rho);

  auto check_rho_nonnegative
      = check_cl(function, "rho", rho_val, "nonnegative");
  auto rho_nonnegative = rho_val >= 0.0;
  auto check_gamma_nonnegative
      = check_cl(function, "Gamma", Gamma_val, "nonnegative");
  auto gamma_nonnegative = gamma_val >= 0.0;
  auto gamma_rowwise_sum = rowwise_sum(gamma);
  auto check_gamma_sum_1 = check_cl(function, "Rowwise sum of gamma",
                                    gamma_rowwise_sum, "equal to 1");
  auto gamma_sum_1 = fabs(gamma_rowwise_sum - 1) <= CONSTRAINT_TOLERANCE;
  auto rho_sum = colwise_sum(rho);

  auto omegas = exp(log_omegas);
  auto alphas;
  auto alpha_log_norms;
  auto norm_norm;

  auto inv_sigma = elt_divide(1., sigma_val);
  auto y_scaled = elt_multiply((y_val - mu_val), inv_sigma);
  auto y_scaled_sq = elt_multiply(y_scaled, y_scaled);

  auto logp1 = -0.5 * y_scaled_sq;
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_rho_cl>::value>(
          logp1 - log(sigma_val), logp1));

  auto scaled_diff = elt_multiply(inv_sigma, y_scaled);
  auto sigma_deriv = elt_multiply(inv_sigma, y_scaled_sq) - inv_sigma;

  matrix_cl<double> logp_cl;
  matrix_cl<double> mu_deriv_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> sigma_deriv_cl;

  results(check_y_not_nan, check_mu_finite, check_sigma_positive, logp_cl,
          y_deriv_cl, mu_deriv_cl, sigma_deriv_cl)
      = expressions(y_not_nan, mu_finite, sigma_positive, logp_expr,
                    calc_if<!is_constant<T_omega_cl>::value>(-scaled_diff),
                    calc_if<!is_constant<T_Gamma_cl>::value>(scaled_diff),
                    calc_if<!is_constant<T_rho_cl>::value>(sigma_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }

  operands_and_partials<T_omega_cl, T_Gamma_cl, T_rho_cl> ops_partials(y, mu,
                                                                       sigma);

  if (!is_constant<T_omega_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_Gamma_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(mu_deriv_cl);
  }
  if (!is_constant<T_rho_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(sigma_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
