#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lccdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::exp;
  using std::log;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "normal_lccdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  T_partials_return ccdf_log(0.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(y, mu, sigma);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return mu_dbl = mu_vec.val(n);
    const T_partials_return sigma_dbl = sigma_vec.val(n);

    const T_partials_return scaled_diff
        = (y_dbl - mu_dbl) / (sigma_dbl * SQRT_TWO);

    T_partials_return one_m_erf;
    if (scaled_diff < -37.5 * INV_SQRT_TWO) {
      one_m_erf = 2.0;
    } else if (scaled_diff < -5.0 * INV_SQRT_TWO) {
      one_m_erf = 2.0 - erfc(-scaled_diff);
    } else if (scaled_diff > 8.25 * INV_SQRT_TWO) {
      one_m_erf = 0.0;
    } else {
      one_m_erf = 1.0 - erf(scaled_diff);
    }

    ccdf_log += LOG_HALF + log(one_m_erf);

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      const T_partials_return rep_deriv_div_sigma
          = scaled_diff > 8.25 * INV_SQRT_TWO
                ? INFTY
                : SQRT_TWO_OVER_SQRT_PI * exp(-scaled_diff * scaled_diff)
                      / one_m_erf / sigma_dbl;
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[n] -= rep_deriv_div_sigma;
      }
      if (!is_constant_all<T_loc>::value) {
        partials<1>(ops_partials)[n] += rep_deriv_div_sigma;
      }
      if (!is_constant_all<T_scale>::value) {
        partials<2>(ops_partials)[n]
            += rep_deriv_div_sigma * scaled_diff * SQRT_TWO;
      }
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
