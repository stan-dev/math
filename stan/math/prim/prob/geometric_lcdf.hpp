#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_LCDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log CDF of the geometric distribution with success probability
 * parameter. Given containers of matching sizes, returns the log sum of
 * probabilities.
 *
 * log CDF: log(1 - (1 - theta)^(n + 1))
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if theta is not in [0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_lcdf(const T_n& n,
                                            const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_prob_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_lcdf";

  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  if (size_zero(n, theta)) {
    return 0.0;
  }

  T_n_ref n_ref = n;
  T_prob_ref theta_ref = theta;
  check_bounded(function, "Success probability parameter",
                value_of(theta_ref), 0.0, 1.0);

  scalar_seq_view<T_n_ref> n_vec(n_ref);
  scalar_seq_view<T_prob_ref> theta_vec(theta_ref);
  const size_t max_size_seq_view = max_size(n_ref, theta_ref);

  for (int i = 0; i < stan::math::size(n); i++) {
    if (n_vec.val(i) < 0) {
      return negative_infinity();
    }
  }

  T_partials_return log_cdf(0.0);
  auto ops_partials = make_partials_propagator(theta_ref);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    const auto theta_val = theta_vec.val(i);
    const auto n_val = n_vec.val(i);
    const auto np1 = n_val + 1.0;
    const auto log1m_theta = log1m(theta_val);
    const auto log_ccdf = np1 * log1m_theta;

    log_cdf += log1m_exp(log_ccdf);

    if constexpr (is_autodiff_v<T_prob>) {
      const auto ccdf = stan::math::exp(log_ccdf);
      partials<0>(ops_partials)[i]
          += np1 * ccdf / ((1.0 - theta_val) * (1.0 - ccdf));
    }
  }

  return ops_partials.build(log_cdf);
}

}  // namespace math
}  // namespace stan
#endif
