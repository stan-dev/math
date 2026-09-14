#ifndef STAN_MATH_PRIM_PROB_LOGISTIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGISTIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

// Logistic(y|mu, sigma) [sigma > 0]
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> logistic_cdf(const T_y& y,
                                                       const T_loc& mu,
                                                       const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "logistic_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(y, mu, sigma);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == NEGATIVE_INFTY) {
      return ops_partials.build(0.0);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      continue;
    }

    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return mu_dbl = mu_vec.val(n);
    const T_partials_return sigma_inv_vec = 1.0 / sigma_vec.val(n);
    const T_partials_return scaled_diff = (y_dbl - mu_dbl) * sigma_inv_vec;
    const T_partials_return Pn = inv_logit(scaled_diff);

    P *= Pn;

    // The partials accumulate d/d. log(Pn); they are rescaled by the product
    // P below. Writing the log-scale derivative as inv_logit(-scaled_diff)
    // avoids the pdf / Pn quotient, which is 0 / 0 once Pn underflows.
    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
      const T_partials_return deriv = inv_logit(-scaled_diff) * sigma_inv_vec;
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += deriv;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= deriv;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n] -= scaled_diff * deriv;
      }
    }
  }

  if constexpr (is_autodiff_v<T_y>) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_loc>) {
    for (size_t n = 0; n < stan::math::size(mu); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  if constexpr (is_autodiff_v<T_scale>) {
    for (size_t n = 0; n < stan::math::size(sigma); ++n) {
      partials<2>(ops_partials)[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
