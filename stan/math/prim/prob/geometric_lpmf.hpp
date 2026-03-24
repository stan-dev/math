#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the geometric distribution. If containers
 * of matching sizes are supplied, returns the log sum of probabilities.
 *
 * The geometric distribution with success probability \f$\theta\f$ has PMF
 * \f[
 *   P(n \mid \theta) = \theta (1 - \theta)^{n}
 * \f]
 * where \f$n \in \{0, 1, 2, \ldots\}\f$ and
 *       \f$\theta \in (0, 1]\f$.
 *
 * This is the number-of-failures parameterization, consistent with
 * the geometric distribution being a special case of the negative
 * binomial distribution with r = 1.
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not in [0, 1]
 * @throw std::domain_error if n is negative
 * @throw std::invalid_argument if container sizes mismatch
 */
template <bool propto, typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n,
                                            const T_prob& theta) {
  using std::log;
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_prob_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lpmf";
  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  if (size_zero(n, theta)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_prob_ref theta_ref = theta;
  check_nonnegative(function, "Outcome variable", n_ref);
  check_bounded(function, "Success probability parameter",
                value_of(theta_ref), 0.0, 1.0);

  if constexpr (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(theta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_prob_ref> theta_vec(theta_ref);
  const size_t max_size_seq_view = max_size(n_ref, theta_ref);
  T_partials_return logp(0.0);

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const auto theta_val = theta_vec.val(i);
    const auto n_val = n_vec.val(i);

    // When theta == 1.0, P(n=0) = 1, P(n>0) = 0
    if (theta_val == 1.0) {
      if (n_val > 0) {
        return negative_infinity();
      }
      // logp += 0 for n=0, theta=1
      continue;
    }

    logp += log(theta_val) + n_val * log1m(theta_val);

    if constexpr (is_autodiff_v<T_prob>) {
      partials<0>(ops_partials)[i]
          += 1.0 / theta_val - n_val / (1.0 - theta_val);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_lpmf(const T_n& n,
                                            const T_prob& theta) {
  return geometric_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
