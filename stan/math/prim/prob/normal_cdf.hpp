#ifndef STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the normal cumulative distribution function for the given
 * variate, location, and scale.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * @tparam T_y type of y
 * @tparam T_loc type of mean parameter
 * @tparam T_scale type of standard deviation parameter
 * @param y A scalar variate.
 * @param mu The location of the normal distribution.
 * @param sigma The scale of the normal distribution
 * @return The unit normal cdf evaluated at the specified arguments.
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_cdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::exp;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "normal_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 1.0;
  }

  T_partials_return cdf(1.0);
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

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
    T_partials_return cdf_n;
    if (scaled_diff < -37.5 * INV_SQRT_TWO) {
      cdf_n = 0.0;
    } else if (scaled_diff < -5.0 * INV_SQRT_TWO) {
      cdf_n = 0.5 * erfc(-scaled_diff);
    } else if (scaled_diff > 8.25 * INV_SQRT_TWO) {
      cdf_n = 1;
    } else {
      cdf_n = 0.5 * (1.0 + erf(scaled_diff));
    }

    cdf *= cdf_n;

    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      const T_partials_return rep_deriv
          = (scaled_diff < -37.5 * INV_SQRT_TWO)
                ? 0.0
                : INV_SQRT_TWO_PI * exp(-scaled_diff * scaled_diff)
                      / (cdf_n * sigma_dbl);
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n] += rep_deriv;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_[n] -= rep_deriv;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge3_.partials_[n] -= rep_deriv * scaled_diff * SQRT_TWO;
      }
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= cdf;
    }
  }
  if (!is_constant_all<T_loc>::value) {
    for (size_t n = 0; n < stan::math::size(mu); ++n) {
      ops_partials.edge2_.partials_[n] *= cdf;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < stan::math::size(sigma); ++n) {
      ops_partials.edge3_.partials_[n] *= cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
