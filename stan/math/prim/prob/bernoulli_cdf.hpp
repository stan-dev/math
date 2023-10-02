#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_CDF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the Bernoulli distribution. If containers are
 * supplied, returns the product of the probabilities.
 *
 * @tparam T_n type of integer parameter
 * @tparam T_prob type of chance of success parameter
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
return_type_t<T_prob> bernoulli_cdf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  static const char* function = "bernoulli_cdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_theta_ref theta_ref = theta;
  const auto& n_arr = as_value_column_array_or_scalar(n);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);
  check_bounded(function, "Probability parameter", theta_arr, 0.0, 1.0);

  if (size_zero(n, theta)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  if (any(n_arr < 0)) {
    return ops_partials.build(0.0);
  }
  const auto& log1m_theta = select(theta_arr == 1, 0.0, log1m(theta_arr));
  const auto& P1 = select(n_arr == 0, log1m_theta, 0.0);

  T_partials_return P = sum(P1);

  if (!is_constant_all<T_prob>::value) {
    partials<0>(ops_partials) = select(n_arr == 0, -exp(P - P1), 0.0);
  }
  return ops_partials.build(exp(P));
}

}  // namespace math
}  // namespace stan
#endif
