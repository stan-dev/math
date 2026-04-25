#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the geometric distribution. Given containers of
 * matching sizes, returns the log of the product of complementary
 * probabilities.
 *
 * log P(N > n | theta) = log((1 - theta)^(n + 1)) = (n + 1) * log1m(theta).
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log complementary probability or log product of complements
 * @throw std::domain_error if theta is not in (0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_lccdf(const T_n& n,
                                             const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lccdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_theta_ref theta_ref = theta;
  const auto& n_arr = as_value_column_array_or_scalar(n);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);
  check_positive_finite(function, "Probability parameter", theta_arr);
  check_less_or_equal(function, "Probability parameter", theta_arr, 1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // log P(N > n) = 0 (i.e. P = 1) when n < 0, matching the existing
  // implementation that short-circuits on the first negative element.
  if (any(n_arr < 0)) {
    return ops_partials.build(0.0);
  }

  // n at INT_MAX: P(N > n) underflows to 0, lccdf = -inf.
  // (The autodiff test framework probes the upper bound at INT_MAX,
  // mirroring the early return used in neg_binomial_lccdf.)
  if (any(n_arr == std::numeric_limits<int>::max())) {
    return ops_partials.build(NEGATIVE_INFTY);
  }

  // theta = 1 means certain success, so P(N > n) = 0 for n >= 0 and the
  // log is -inf. The partials path divides by (theta - 1) = 0, so we
  // short-circuit.
  if (any(theta_arr == 1.0)) {
    return ops_partials.build(NEGATIVE_INFTY);
  }

  // log P(N > n) = (n + 1) * log1m(theta)
  const auto& log1m_theta = log1m(theta_arr);
  T_partials_return logP = sum((n_arr + 1.0) * log1m_theta);

  if constexpr (is_autodiff_v<T_prob>) {
    // d/dtheta (n + 1) * log1m(theta) = -(n + 1) / (1 - theta)
    //                                 = (n + 1) / (theta - 1)
    // theta = 1 case was filtered above so theta - 1 != 0 here.
    if constexpr (is_stan_scalar_v<T_prob>) {
      partials<0>(ops_partials) = sum((n_arr + 1.0) * inv(theta_arr - 1.0));
    } else {
      partials<0>(ops_partials) = (n_arr + 1.0) * inv(theta_arr - 1.0);
    }
  }

  return ops_partials.build(logP);
}

}  // namespace math
}  // namespace stan
#endif
