#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the logit-parametrized Bernoulli distribution. If
 * containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n type of integer parameter
 * @tparam T_prob type of chance of success parameter
 * @param n integer parameter
 * @param theta logit-transformed chance of success parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is infinite.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n, typename T_prob>
return_type_t<T_prob> bernoulli_logit_lpmf(const T_n& n, const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using std::exp;
  static const char* function = "bernoulli_logit_lpmf";
  check_bounded(function, "n", n, 0, 1);
  check_not_nan(function, "Logit transformed probability parameter", theta);
  check_consistent_sizes(function, "Random variable", n,
                         "Probability parameter", theta);

  if (size_zero(n, theta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_prob>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_prob> ops_partials(theta);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_prob> theta_vec(theta);
  size_t N = max_size(n, theta);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return theta_dbl = value_of(theta_vec[n]);

    const int sign = 2 * n_vec[n] - 1;
    const T_partials_return ntheta = sign * theta_dbl;
    const T_partials_return exp_m_ntheta = exp(-ntheta);

    // Handle extreme values gracefully using Taylor approximations.
    static const double cutoff = 20.0;
    if (ntheta > cutoff) {
      logp -= exp_m_ntheta;
    } else if (ntheta < -cutoff) {
      logp += ntheta;
    } else {
      logp -= log1p(exp_m_ntheta);
    }

    if (!is_constant_all<T_prob>::value) {
      if (ntheta > cutoff) {
        ops_partials.edge1_.partials_[n] -= exp_m_ntheta;
      } else if (ntheta < -cutoff) {
        ops_partials.edge1_.partials_[n] += sign;
      } else {
        ops_partials.edge1_.partials_[n]
            += sign * exp_m_ntheta / (exp_m_ntheta + 1);
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_prob>
inline return_type_t<T_prob> bernoulli_logit_lpmf(const T_n& n,
                                                  const T_prob& theta) {
  return bernoulli_logit_lpmf<false>(n, theta);
}

}  // namespace math
}  // namespace stan
#endif
