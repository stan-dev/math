#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_LPDF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

/**
 * From this value of nu on, student_t_lpdf forms its nu terms with the
 * asymptotic series below. Below it, the differences of two lgamma and of
 * two digamma values lose about eps * lgamma(nu / 2), at most about 1e-14.
 * For large nu that loss grows like eps * nu log(nu), and the result has no
 * correct digits from about nu = 1e16. At nu = 40 the
 * first omitted series term is below 1e-17 of the result.
 */
constexpr double student_t_series_min_nu = 40.0;

/**
 * Return lgamma((nu + 1) / 2) - lgamma(nu / 2) - log(nu) / 2, the part of
 * the Student-t normalizing constant that depends on nu. For nu >= 40,
 * with h = nu / 2, the asymptotic series
 *
 *   lgamma(h + 1/2) - lgamma(h) = log(h) / 2 - 1 / (8 h) + 1 / (192 h^3)
 *       - 1 / (640 h^5) + 17 / (14336 h^7) - 31 / (18432 h^9) + O(h^-11)
 *
 * is used, in which log(h) / 2 - log(nu) / 2 = -log(2) / 2.
 */
struct student_t_nu_constant_fun {
  template <typename T>
  static inline T fun(const T& nu) {
    if (value_of_rec(nu) < student_t_series_min_nu) {
      const T half_nu = 0.5 * nu;
      return lgamma(half_nu + 0.5) - lgamma(half_nu) - 0.5 * log(nu);
    }
    const T z = 2.0 / nu;
    const T z2 = square(z);
    return -0.5 * LOG_TWO
           + z
                 * (-0.125
                    + z2
                          * (1.0 / 192.0
                             + z2
                                   * (-1.0 / 640.0
                                      + z2
                                            * (17.0 / 14336.0
                                               - z2 * (31.0 / 18432.0)))));
  }
};

/**
 * Return digamma((nu + 1) / 2) - digamma(nu / 2). For nu >= 40, with
 * h = nu / 2, the derivative in h of the series above is used:
 *
 *   1 / (2 h) + 1 / (8 h^2) - 1 / (64 h^4) + 1 / (128 h^6)
 *       - 119 / (14336 h^8) + 279 / (18432 h^10) + O(h^-12).
 */
struct student_t_nu_digamma_fun {
  template <typename T>
  static inline T fun(const T& nu) {
    if (value_of_rec(nu) < student_t_series_min_nu) {
      const T half_nu = 0.5 * nu;
      return digamma(half_nu + 0.5) - digamma(half_nu);
    }
    const T z = 2.0 / nu;
    const T z2 = square(z);
    return z
           * (0.5
              + z
                    * (0.125
                       + z2
                             * (-1.0 / 64.0
                                + z2
                                      * (1.0 / 128.0
                                         + z2
                                               * (-119.0 / 14336.0
                                                  + z2 * (279.0 / 18432.0))))));
  }
};

}  // namespace internal

/** \ingroup prob_dists
 * The log of the Student-t density for the given y, nu, mean, and
 * scale parameter.  The scale parameter must be greater
 * than 0.
 *
 * \f{eqnarray*}{
 y &\sim& t_{\nu} (\mu, \sigma^2) \\
 \log (p (y \, |\, \nu, \mu, \sigma) ) &=& \log \left( \frac{\Gamma((\nu + 1)
 /2)}
 {\Gamma(\nu/2)\sqrt{\nu \pi} \sigma} \left( 1 + \frac{1}{\nu} (\frac{y -
 \mu}{\sigma})^2 \right)^{-(\nu + 1)/2} \right) \\
 &=& \log( \Gamma( (\nu+1)/2 )) - \log (\Gamma (\nu/2) - \frac{1}{2} \log(\nu
 \pi) - \log(\sigma)
 -\frac{\nu + 1}{2} \log (1 + \frac{1}{\nu} (\frac{y - \mu}{\sigma})^2)
 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_dof type of degrees of freedom
 * @tparam T_loc type of location
 * @tparam T_scale type of scale
 *
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @param mu The mean of the Student-t distribution.
 * @param sigma The scale parameter of the Student-t distribution.
 * @return The log of the Student-t density at y.
 * @throw std::domain_error if sigma is not greater than 0.
 * @throw std::domain_error if nu is not greater than 0.
 */
template <bool propto, typename T_y, typename T_dof, typename T_loc,
          typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_dof, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_dof, T_loc, T_scale> student_t_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_nu_ref = ref_type_if_not_constant_t<T_dof>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "student_t_lpdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu,
                         "Location parameter", mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_nu_ref nu_ref = nu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) nu_val = to_ref(as_value_column_array_or_scalar(nu_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_not_nan(function, "Random variable", y_val);
  check_positive_finite(function, "Degrees of freedom parameter", nu_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, nu, mu, sigma)) {
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_y, T_dof, T_loc, T_scale>::value) {
    return 0.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, nu_ref, mu_ref, sigma_ref);

  const auto& half_nu
      = to_ref_if<include_summand<propto, T_dof>::value>(0.5 * nu_val);
  const auto& square_y_scaled = square((y_val - mu_val) / sigma_val);
  const auto& square_y_scaled_over_nu
      = to_ref_if<is_any_autodiff_v<T_y, T_dof, T_loc, T_scale>>(square_y_scaled
                                                                 / nu_val);
  const auto& log1p_val
      = to_ref_if<is_autodiff_v<T_dof>>(log1p(square_y_scaled_over_nu));

  size_t N = max_size(y, nu, mu, sigma);
  T_partials_return logp = -sum((half_nu + 0.5) * log1p_val);
  if constexpr (include_summand<propto>::value) {
    logp -= LOG_SQRT_PI * N;
  }
  if constexpr (include_summand<propto, T_dof>::value) {
    // lgamma(nu/2 + 1/2) - lgamma(nu/2) - log(nu)/2, by a series for large
    // nu; see internal::student_t_nu_constant_fun
    logp += sum(apply_scalar_unary<
                internal::student_t_nu_constant_fun,
                std::decay_t<decltype(nu_val)>>::apply(nu_val))
            * N / math::size(nu);
  }
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }

  if constexpr (is_any_autodiff_v<T_y, T_loc>) {
    const auto& square_sigma = square(sigma_val);
    auto deriv_y_mu = to_ref_if<(is_autodiff_v<T_y> && is_autodiff_v<T_loc>)>(
        (nu_val + 1) * (y_val - mu_val)
        / ((1 + square_y_scaled_over_nu) * square_sigma * nu_val));
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = -deriv_y_mu;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<2>(ops_partials) = std::move(deriv_y_mu);
    }
  }
  if constexpr (is_any_autodiff_v<T_dof, T_scale>) {
    const auto& rep_deriv
        = to_ref_if<(is_autodiff_v<T_dof> && is_autodiff_v<T_scale>)>(
            (nu_val + 1) * square_y_scaled_over_nu
                / (1 + square_y_scaled_over_nu)
            - 1);
    if constexpr (is_autodiff_v<T_dof>) {
      // digamma(nu/2 + 1/2) - digamma(nu/2), by a series for large nu;
      // see internal::student_t_nu_digamma_fun
      const auto& digamma_nu_term
          = apply_scalar_unary<internal::student_t_nu_digamma_fun,
                               std::decay_t<decltype(nu_val)>>::apply(nu_val);
      edge<1>(ops_partials).partials_
          = 0.5 * (digamma_nu_term - log1p_val + rep_deriv / nu_val);
    }
    if constexpr (is_autodiff_v<T_scale>) {
      partials<3>(ops_partials) = rep_deriv / sigma_val;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_dof, T_loc, T_scale> student_t_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& sigma) {
  return student_t_lpdf<false>(y, nu, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
