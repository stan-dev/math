#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/select.hpp>
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
  using T_n_ref = ref_type_t<T_n>;
  using T_theta_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lccdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_n_ref n_ref = n;
  T_theta_ref theta_ref = theta;
  const auto& n_arr = as_value_column_array_or_scalar(n_ref);
  const auto& theta_arr = as_value_column_array_or_scalar(theta_ref);
  check_positive_finite(function, "Probability parameter", theta_arr);
  check_less_or_equal(function, "Probability parameter", theta_arr, 1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  // n at INT_MAX: P(N > n) underflows to 0, lccdf = -inf.
  // (The autodiff test framework probes the upper bound at INT_MAX,
  // mirroring the early return used in neg_binomial_lccdf.)
  if (any(n_arr == std::numeric_limits<int>::max())) {
    return ops_partials.build(NEGATIVE_INFTY);
  }

  // theta = 1 with any n_i >= 0 makes (n_i + 1) * log1m(1) = -inf and
  // the partials path divides by (theta - 1) = 0. Per-element loop here
  // correctly pairs theta_i with n_i (e.g. theta=[0.5, 1.0], n=[2, -1]
  // returns 3*log(0.5), not -inf). For theta_i = 1 paired with n_i < 0
  // the per-element select below contributes 0 to both value and partials.
  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t max_sz = max_size(n_ref, theta_ref);
  for (size_t i = 0; i < max_sz; i++) {
    if (value_of(theta_vec[i]) == 1.0 && n_vec[i] >= 0) {
      return ops_partials.build(NEGATIVE_INFTY);
    }
  }

  // Per-element: n_i < 0 contributes log(1) = 0 (since P(N > n_i) = 1
  // for any n_i below the support), so they are no-ops in the sum.
  // n_i >= 0 contributes (n_i + 1) * log1m(theta).
  const auto& log1m_theta = log1m(theta_arr);
  const auto& term
      = select(n_arr < 0, T_partials_return(0), (n_arr + 1.0) * log1m_theta);
  T_partials_return logP = sum(term);

  if constexpr (is_autodiff_v<T_prob>) {
    // d/dtheta of 0 = 0 for n_i < 0, and (n + 1) / (theta - 1) for n_i >= 0.
    // theta = 1 was filtered above, so theta - 1 != 0 here.
    const auto& dterm = select(n_arr < 0, T_partials_return(0),
                               (n_arr + 1.0) * inv(theta_arr - 1.0));
    if constexpr (is_stan_scalar_v<T_prob>) {
      partials<0>(ops_partials) = sum(dterm);
    } else {
      partials<0>(ops_partials) = dterm;
    }
  }

  return ops_partials.build(logP);
}

}  // namespace math
}  // namespace stan
#endif
