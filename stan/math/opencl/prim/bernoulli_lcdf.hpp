#ifndef STAN_MATH_OPENCL_PRIM_BERNOULLI_LCDF_HPP
#define STAN_MATH_OPENCL_PRIM_BERNOULLI_LCDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CDF of the Bernoulli distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n_cl type of integer parameter
 * @tparam T_prob_cl type of chance of success parameter
 * @param n integer parameter
 * @param theta logit-transformed chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <
    typename T_n_cl, typename T_prob_cl,
    require_all_prim_or_rev_kernel_expression_t<T_n_cl, T_prob_cl>* = nullptr,
    require_any_not_stan_scalar_t<T_n_cl, T_prob_cl>* = nullptr>
return_type_t<T_prob_cl> bernoulli_lcdf(const T_n_cl& n,
                                        const T_prob_cl& theta) {
  static const char* function = "bernoulli_lcdf(OpenCL)";
  using T_partials_return = partials_return_t<T_prob_cl>;
  using std::isnan;
  constexpr bool is_n_vector = !is_stan_scalar<T_n_cl>::value;
  constexpr bool is_theta_vector = !is_stan_scalar<T_prob_cl>::value;

  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  const size_t N = is_n_vector ? math::size(n) : math::size(theta);
  if (N == 0) {
    return 0.0;
  }

  const auto& theta_col = as_column_vector_or_scalar(theta);
  const auto& theta_val = value_of(theta_col);

  auto check_theta_bounded = check_cl(function, "Probability parameter",
                                      theta_val, "in the interval [0, 1]");
  auto theta_bounded_expr = 0.0 <= theta_val && theta_val <= 1.0;

  auto any_n_negative = colwise_max(cast<char>(n < 0));
  auto Pi = 1.0 - theta_val;
  auto cond = n >= 1;
  auto P_expr = colwise_sum(select(cond, 0.0, log(Pi)));
  auto deriv = select(cond, 0.0, elt_divide(-1.0, Pi));

  matrix_cl<char> any_n_negative_cl;
  matrix_cl<double> P_cl;
  matrix_cl<double> deriv_cl;

  results(check_theta_bounded, any_n_negative_cl, P_cl, deriv_cl)
      = expressions(theta_bounded_expr, any_n_negative, P_expr,
                    calc_if<(!is_constant_all<T_prob_cl>::value)>(deriv));

  if (from_matrix_cl(any_n_negative_cl).maxCoeff()) {
    return NEGATIVE_INFTY;
  }

  T_partials_return P = from_matrix_cl(P_cl).sum();
  operands_and_partials<decltype(theta_col)> ops_partials(theta_col);

  if (!is_constant_all<T_prob_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(deriv_cl);
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
#endif
