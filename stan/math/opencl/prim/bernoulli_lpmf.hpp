#ifndef STAN_MATH_OPENCL_PRIM_BERNOULLI_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_BERNOULLI_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the log PMF of the Bernoulli distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n_cl type of integer parameters
 * @tparam T_prob_cl type of chance of success parameters
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <
    bool propto, typename T_n_cl, typename T_prob_cl,
    require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_prob_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_n_cl, T_prob_cl>* = nullptr>
return_type_t<T_prob_cl> bernoulli_lpmf(const T_n_cl& n,
                                        const T_prob_cl& theta) {
  static const char* function = "bernoulli_lpmf(OpenCL)";
  using T_partials_return = partials_return_t<T_prob_cl>;
  constexpr bool is_n_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_theta_vector = !is_stan_scalar<T_prob_cl>::value;

  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  const size_t N = is_n_vector ? math::size(n) : math::size(theta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob_cl>::value) {
    return 0.0;
  }

  const auto& theta_col = as_column_vector_or_scalar(theta);
  const auto& theta_val = value_of(theta_col);

  T_partials_return logp(0.0);
  operands_and_partials<decltype(theta_col)> ops_partials(theta_col);

  auto check_n_bounded = check_cl(function, "n", n, "in the interval [0, 1]");
  auto n_bounded_expr = 0 <= n && n <= 1;

  if (is_theta_vector) {
    auto logp_expr
        = colwise_sum(select(n == 1, log(theta_val), log1p(-theta_val)));
    auto deriv_expr = inv(theta_val + select(n == 1, 0, -1));

    auto check_theta_bounded = check_cl(function, "Probability parameter",
                                        theta_val, "in the interval [0, 1]");
    auto theta_bounded_expr = 0 <= theta_val && theta_val <= 1;

    matrix_cl<double> logp_cl;
    matrix_cl<double> deriv_cl;

    results(logp_cl, deriv_cl, check_n_bounded, check_theta_bounded)
        = expressions(logp_expr,
                      calc_if<!is_constant_all<T_prob_cl>::value>(deriv_expr),
                      n_bounded_expr, theta_bounded_expr);

    logp = sum(from_matrix_cl(logp_cl));

    if (!is_constant_all<T_prob_cl>::value) {
      ops_partials.edge1_.partials_ = deriv_cl;
    }
  } else {
    auto n_sum_expr = rowwise_sum(forward_as<matrix_cl<int>>(n));

    matrix_cl<int> n_sum_cl;

    results(n_sum_cl, check_n_bounded)
        = expressions(n_sum_expr, n_bounded_expr);

    size_t n_sum = sum(from_matrix_cl(n_sum_cl));
    double theta_val_scal = forward_as<double>(theta_val);
    if (n_sum == N) {
      logp = N * log(theta_val_scal);
    } else if (n_sum == 0) {
      logp = N * log1m(theta_val_scal);
    } else {
      logp = n_sum * log(theta_val_scal) + (N - n_sum) * log1m(theta_val_scal);
    }
    if (!is_constant_all<T_prob_cl>::value) {
      double& edge1_partial = forward_as<internal::broadcast_array<double>>(
          ops_partials.edge1_.partials_)[0];
      if (n_sum == N) {
        edge1_partial += N / theta_val_scal;
      } else if (n_sum == 0) {
        edge1_partial += N / (theta_val_scal - 1);
      } else {
        edge1_partial
            += n_sum / theta_val_scal + (N - n_sum) / (theta_val_scal - 1);
      }
    }
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
