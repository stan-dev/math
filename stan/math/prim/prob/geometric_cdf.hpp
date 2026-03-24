#ifndef STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_GEOMETRIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the geometric distribution with success probability
 * parameter. Given containers of matching sizes, returns the product of
 * probabilities.
 *
 * CDF: F(n | theta) = 1 - (1 - theta)^(n + 1)
 *
 * @tparam T_n type of outcome variable
 * @tparam T_prob type of success probability parameter
 *
 * @param n outcome variable (number of failures before first success)
 * @param theta success probability parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if theta is not in [0, 1]
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_prob>
inline return_type_t<T_prob> geometric_cdf(const T_n& n,
                                           const T_prob& theta) {
  using T_partials_return = partials_return_t<T_n, T_prob>;
  using T_n_ref = ref_type_t<T_n>;
  using T_prob_ref = ref_type_t<T_prob>;
  static constexpr const char* function = "geometric_cdf";

  check_consistent_sizes(function, "Outcome variable", n,
                         "Success probability parameter", theta);
  if (size_zero(n, theta)) {
    return 1.0;
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
      return 0.0;
    }
  }

  T_partials_return cdf(1.0);
  auto ops_partials = make_partials_propagator(theta_ref);
  for (size_t i = 0; i < max_size_seq_view; i++) {
    const auto theta_val = theta_vec.val(i);
    const auto n_val = n_vec.val(i);
    const auto np1 = n_val + 1.0;
    const auto log1m_theta = log1m(theta_val);
    const auto ccdf_i = stan::math::exp(np1 * log1m_theta);
    const auto cdf_i = 1.0 - ccdf_i;

    cdf *= cdf_i;

    if constexpr (is_autodiff_v<T_prob>) {
      if (cdf_i > 0.0 && theta_val > 0.0) {
        partials<0>(ops_partials)[i]
            += np1 * ccdf_i / ((1.0 - theta_val) * cdf_i);
      }
    }
  }

  if constexpr (is_autodiff_v<T_prob>) {
    for (size_t i = 0; i < stan::math::size(theta); ++i) {
      partials<0>(ops_partials)[i] *= cdf;
    }
  }

  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
