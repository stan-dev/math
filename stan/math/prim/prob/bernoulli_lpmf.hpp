#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the Bernoulli distribution. If containers are
 * supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n type of integer parameters
 * @tparam T_prob type of chance of success parameters
 * @param n integer parameter
 * @param theta chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n, typename T_prob,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_prob>* = nullptr>
return_type_t<T_prob> bernoulli_lpmf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_theta_ref = ref_type_t<T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using std::log;
  static constexpr const char* function = "bernoulli_lpmf";
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);
  const T_n_ref n_ref = to_ref(n);
  const T_theta_ref theta_ref = to_ref(theta);
  check_bounded(function, "n", n_ref, 0, 1);
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  if (size_zero(n, theta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  auto ops_partials = make_partials_propagator(theta_ref);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t N = max_size(n, theta);

  if (math::size(theta) == 1) {
    size_t sum = 0;
    for (size_t n = 0; n < N; n++) {
      sum += n_vec.val(n);
    }
    const T_partials_return theta_dbl = theta_vec.val(0);
    // avoid nans when sum == N or sum == 0
    if (sum == N) {
      logp += N * log(theta_dbl);
      if (!is_constant_all<T_prob>::value) {
        partials<0>(ops_partials)[0] += N / theta_dbl;
      }
    } else if (sum == 0) {
      logp += N * log1m(theta_dbl);
      if (!is_constant_all<T_prob>::value) {
        partials<0>(ops_partials)[0] += N / (theta_dbl - 1);
      }
    } else {
      const T_partials_return log_theta = log(theta_dbl);
      const T_partials_return log1m_theta = log1m(theta_dbl);

      logp += sum * log_theta;
      logp += (N - sum) * log1m_theta;

      if (!is_constant_all<T_prob>::value) {
        partials<0>(ops_partials)[0] += sum * inv(theta_dbl);
        partials<0>(ops_partials)[0] += (N - sum) * inv(theta_dbl - 1);
      }
    }
  } else {
    for (size_t n = 0; n < N; n++) {
      const int n_int = n_vec.val(n);
      const T_partials_return theta_dbl = theta_vec.val(n);

      if (n_int == 1) {
        logp += log(theta_dbl);
      } else {
        logp += log1m(theta_dbl);
      }

      if (!is_constant_all<T_prob>::value) {
        if (n_int == 1) {
          partials<0>(ops_partials)[n] += inv(theta_dbl);
        } else {
          partials<0>(ops_partials)[n] += inv(theta_dbl - 1);
        }
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_prob>
inline return_type_t<T_prob> bernoulli_lpmf(const T_y& n, const T_prob& theta) {
  return bernoulli_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
