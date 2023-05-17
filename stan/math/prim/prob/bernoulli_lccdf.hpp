#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CCDF of the Bernoulli distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n type of integer parameter
 * @tparam T_prob type of chance of success parameter
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
return_type_t<T_prob> bernoulli_lccdf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  using std::log;
  static const char* function = "bernoulli_lccdf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  T_theta_ref theta_ref = theta;
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(theta_ref);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t max_size_seq_view = max_size(n, theta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(n); i++) {
    const double n_dbl = n_vec.val(i);
    if (n_dbl < 0) {
      return ops_partials.build(0.0);
    }
    if (n_dbl >= 1) {
      return ops_partials.build(NEGATIVE_INFTY);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    const T_partials_return Pi = theta_vec.val(i);

    P += log(Pi);

    if (!is_constant_all<T_prob>::value) {
      partials<0>(ops_partials)[i] += inv(Pi);
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
