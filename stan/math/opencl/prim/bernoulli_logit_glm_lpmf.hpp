#ifndef STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BERNOULLI_LOGIT_GLM_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

#include <cmath>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Bernoulli distribution and logit link function.
 * This is an overload of the GLM in prim/prob/bernoulli_logit_glm_lpmf.hpp
 * that is implemented in OpenCL.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y_cl binary scalar or vector parameter on OpenCL device. If it is a
 * scalar it will be broadcast - used for all instances.
 * @param x_cl design matrix on OpenCL device. This overload does not support
 * broadcasting of a row vector x!
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is not binary.
 * @throw std::invalid_argument if container sizes mismatch.
 */

template <bool propto, typename T_alpha, typename T_beta>
return_type_t<T_alpha, T_beta> bernoulli_logit_glm_lpmf(
    const matrix_cl<int>& y_cl, const matrix_cl<double>& x_cl,
    const T_alpha& alpha, const T_beta& beta) {
  static const char* function = "bernoulli_logit_glm_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_alpha, T_beta>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;
  using Eigen::Dynamic;
  using Eigen::Matrix;

  constexpr int is_alpha_vector = is_vector<T_alpha>::value;

  const size_t N = x_cl.rows();
  const size_t M = x_cl.cols();

  if (y_cl.size() != 1) {
    check_size_match(function, "Rows of ", "x_cl", N, "rows of ", "y_cl",
                     y_cl.rows());
  }
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_alpha_vector) {
    check_size_match(function, "Rows of ", "x_cl", N, "size of ", "alpha",
                     stan::math::size(alpha));
  }

  if (N == 0) {
    return 0;
  }
  if (!include_summand<propto, T_alpha, T_beta>::value) {
    return 0;
  }

  T_beta_ref beta_ref = beta;
  T_alpha_ref alpha_ref = alpha;

  const auto& beta_val = value_of_rec(beta_ref);
  const auto& alpha_val = value_of_rec(alpha_ref);

  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);

  matrix_cl<double> beta_cl(beta_val_vec);
  matrix_cl<double> alpha_cl(alpha_val_vec);

  auto ytheta_expr = matrix_vector_multiply(x_cl, beta_cl)
                     + broadcast<!is_alpha_vector, false>(alpha_cl);
  auto y_bc_expr = colwise_optional_broadcast(y_cl);
  auto signs_expr = 2 * y_bc_expr - 1;
  auto ytheta_signs_expr = elt_multiply(ytheta_expr, signs_expr);
  auto exp_m_ytheta_expr = exp(-ytheta_signs_expr);
  const double cutoff = 20.0;
  auto high_bound_expr = ytheta_signs_expr > cutoff;
  auto low_bound_expr = ytheta_signs_expr < -cutoff;
  auto err_cond_expr = y_bc_expr < 0 || y_bc_expr > 1;
  auto logp_expr
      = colwise_sum(select(err_cond_expr, NOT_A_NUMBER,
                           select(high_bound_expr, -exp_m_ytheta_expr,
                                  select(low_bound_expr, ytheta_signs_expr,
                                         -log1p(exp_m_ytheta_expr)))));
  auto theta_derivative_expr
      = select(high_bound_expr, -exp_m_ytheta_expr,
               select(low_bound_expr, signs_expr,
                      elt_divide(elt_multiply(signs_expr, exp_m_ytheta_expr),
                                 (exp_m_ytheta_expr + 1))));

  const int wgs = logp_expr.rows();
  matrix_cl<double> logp_cl(wgs, 1);
  constexpr bool need_theta_derivative
      = !is_constant_all<T_beta, T_alpha>::value;
  matrix_cl<double> theta_derivative_cl(need_theta_derivative ? N : 0, 1);
  constexpr bool need_theta_derivative_sum
      = need_theta_derivative && !is_alpha_vector;
  matrix_cl<double> theta_derivative_sum_cl(need_theta_derivative_sum ? wgs : 0,
                                            1);

  results(logp_cl, theta_derivative_cl, theta_derivative_sum_cl) = expressions(
      logp_expr, calc_if<need_theta_derivative>(theta_derivative_expr),
      calc_if<need_theta_derivative_sum>(colwise_sum(theta_derivative_expr)));

  T_partials_return logp = sum(from_matrix_cl<Eigen::Dynamic, 1>(logp_cl));
  if (!std::isfinite(logp)) {
    check_bounded(function, "Vector of dependent variables",
                  from_matrix_cl(y_cl), 0, 1);
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables",
                 from_matrix_cl(x_cl));
  }

  operands_and_partials<T_alpha_ref, T_beta_ref> ops_partials(alpha_ref,
                                                              beta_ref);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_alpha>::value) {
    if (is_alpha_vector) {
      ops_partials.edge1_.partials_
          = std::move(from_matrix_cl<Dynamic, 1>(theta_derivative_cl));
    } else {
      ops_partials.edge1_.partials_[0]
          = sum(from_matrix_cl(theta_derivative_sum_cl));
    }
  }
  if (!is_constant_all<T_beta>::value) {
    matrix_cl<double> theta_derivative_transpose_cl(
        theta_derivative_cl.buffer(), 1,
        theta_derivative_cl
            .rows());  // transposition of a vector can be done without copying
    ops_partials.edge2_.partials_
        = from_matrix_cl<1, Dynamic>(theta_derivative_transpose_cl * x_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan

#endif
#endif
