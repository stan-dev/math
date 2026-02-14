#ifndef STAN_MATH_PRIM_PROB_YULE_SIMON_LCDF_HPP
#define STAN_MATH_PRIM_PROB_YULE_SIMON_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CDF of the Yule-Simon distribution with shape parameter.
 * Given containers of matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_n type of outcome parameter
 * @tparam T_alpha type of shape parameter
 *
 * @param n outcome variable
 * @param alpha shape parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if alpha fails to be positive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_alpha>
inline return_type_t<T_alpha> yule_simon_lcdf(const T_n& n,
                                              const T_alpha& alpha) {
  using T_partials_return = partials_return_t<T_n, T_alpha>;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  static constexpr const char* function = "yule_simon_lcdf";

  check_consistent_sizes(function, "Outcome variable", n, "Shape parameter",
                         alpha);
  if (size_zero(n, alpha)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  check_positive_finite(function, "Shape parameter", alpha_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  const size_t max_size_seq_view = max_size(n_ref, alpha_ref);

  // Explicit return for invalid or extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (int i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 1.0) {
      return negative_infinity();
    }
    if (n_vec.val(i) == std::numeric_limits<int>::max()) {
      return 0.0;
    }
  }

  T_partials_return log_cdf(0.0);
  auto ops_partials = make_partials_propagator(alpha_ref);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    auto np1 = n_vec.val(i) + 1.0;
    auto ap1 = alpha_vec.val(i) + 1.0;
    auto nap1 = n_vec.val(i) + ap1;

    auto log_ccdf = lgamma(ap1) + lgamma(np1) - lgamma(nap1);
    log_cdf += log1m_exp(log_ccdf);

    if constexpr (is_autodiff_v<T_alpha>) {
      auto ccdf = stan::math::exp(log_ccdf);
      auto chain_rule_term = -ccdf / (1.0 - ccdf);
      partials<0>(ops_partials)[i]
          += (digamma(ap1) - digamma(nap1)) * chain_rule_term;
    }
  }

  return ops_partials.build(log_cdf);
}

}  // namespace math
}  // namespace stan
#endif
