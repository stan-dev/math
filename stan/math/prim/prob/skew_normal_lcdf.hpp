#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {
constexpr char skew_normal_lcdf_func[] = "skew_normal_lcdf";

/**
 * Log of the standard skew normal cdf F(z; alpha) and its derivatives with
 * respect to z and alpha.
 *
 * @tparam T type of the log cdf and of its derivatives
 */
template <typename T>
struct skew_normal_log_cdf_terms {
  T log_cdf; /**< log F(z; alpha) */
  T d_z;     /**< derivative of log F with respect to z */
  T d_alpha; /**< derivative of log F with respect to alpha */
};

// Nodes and weights of the 32-point Gauss-Laguerre rule, without the 9
// largest nodes, whose weights are below 1e-20.
// clang-format off
static constexpr double skew_normal_gl_nodes[23] = {
    0.04448936583326702, 0.23452610951961853, 0.5768846293018864,
    1.0724487538178176, 1.7224087764446454, 2.5283367064257947,
    3.4922132730219944, 4.616456769749767, 5.903958504174244,
    7.358126733186241, 8.982940924212595, 10.783018632539973,
    12.763697986742725, 14.931139755522556, 17.292454336715316,
    19.855860940336054, 22.630889013196775, 25.628636022459247,
    28.862101816323474, 32.346629153964734, 36.10049480575197,
    40.14571977153944, 44.509207995754934};
static constexpr double skew_normal_gl_weights[23] = {
    0.10921834195238497, 0.21044310793881324, 0.235213229669848,
    0.19590333597288104, 0.12998378628607177, 0.07057862386571744,
    0.03176091250917507, 0.011918214834838558, 0.0037388162946115247,
    0.0009808033066149551, 0.0002148649188013642, 3.920341967987947e-05,
    5.9345416128686326e-06, 7.416404578667552e-07, 7.604567879120781e-08,
    6.350602226625806e-09, 4.281382971040929e-10, 2.3058994918913362e-11,
    9.799379288727094e-13, 3.2378016577292665e-14, 8.171823443420719e-16,
    1.5421338333938235e-17, 2.1197922901636187e-19};
// clang-format on

/**
 * Log of the standard skew normal cdf in the short tail: z < 0, alpha > 0
 * and alpha |z| >= 2.5. There Phi(z) - 2 T(z, alpha) loses about
 * (alpha z)^2 / (2 log 10) digits.
 *
 * With c = z^2 / 2, d1 = (alpha z)^2 / 2 and d2 = c (1 + alpha^2),
 *
 *   F = exp(-c (1 + alpha^2)) J / (2 pi c alpha (1 + alpha^2)),
 *   J = int_0^inf exp(-w) / ((1 + w / d2) sqrt(1 + w / d1)) dw,
 *
 * from F = (1 / pi) int_alpha^inf exp(-c (1 + t^2)) / (1 + t^2) dt and the
 * substitution w = c (t^2 - alpha^2). J lies in (0, 1]. The Gauss-Laguerre
 * rule above gives J with a relative error below 8e-15 for alpha |z| >= 2.5
 * and below 2e-17 for alpha |z| >= 3. The derivatives are
 *
 *   d/dz log F = |z| (1 + alpha^2) M(alpha |z|) / J,
 *   d/dalpha log F = -z^2 alpha / J,
 *
 * with M(x) = sqrt(pi / 2) x erfcx(x / sqrt(2)), which tends to 1.
 *
 * @tparam T type of the arguments
 * @param z standardized variate, z < 0
 * @param alpha shape, alpha > 0
 * @return log F and its derivatives with respect to z and alpha
 */
template <typename T>
inline skew_normal_log_cdf_terms<T> skew_normal_log_cdf_tail(const T& z,
                                                             const T& alpha) {
  const T x = -alpha * z;
  const T z_sq = square(z);
  const T x_sq = square(x);
  const T inv_d1 = 2.0 / x_sq;
  const T inv_d2 = 2.0 / (z_sq + x_sq);
  T j = 0.0;
  for (int i = 0; i < 23; ++i) {
    j += skew_normal_gl_weights[i]
         / ((1.0 + skew_normal_gl_nodes[i] * inv_d2)
            * sqrt(1.0 + skew_normal_gl_nodes[i] * inv_d1));
  }
  const T log_cdf
      = -0.5 * (z_sq + x_sq)
        - (LOG_PI + 2.0 * log(-z) + log(alpha) + log1p(square(alpha))) + log(j);
  // M(x) - 1 is below 1e-16 for x > 1e8; the bound also avoids inf * 0
  const T mills = value_of_rec(x) > 1e8
                      ? T(1.0)
                      : SQRT_PI * INV_SQRT_TWO * x * erfcx(x * INV_SQRT_TWO);
  return {log_cdf, -z * (1.0 + square(alpha)) * mills / j, -z_sq * alpha / j};
}

/**
 * Log of the standard skew normal cdf F(z; alpha) and its derivatives with
 * respect to z and alpha.
 *
 * Three representations, chosen so that no large cancellation remains:
 * the short tail (z < 0, alpha > 0, alpha |z| >= 2.5) uses
 * skew_normal_log_cdf_tail; otherwise alpha > 1 uses
 * F = 2 T(alpha |z|, 1 / alpha) + Phi(alpha z) erf(z / sqrt(2)), from
 * T(h, a) + T(a h, 1 / a) = Phi(h) / 2 + Phi(a h) / 2 - Phi(h) Phi(a h),
 * and alpha <= 1 uses F = Phi(z) - 2 T(z, alpha).
 *
 * @tparam T type of the arguments
 * @param z standardized variate
 * @param alpha shape
 * @return log F and its derivatives with respect to z and alpha
 */
template <typename T>
inline skew_normal_log_cdf_terms<T> skew_normal_log_cdf(const T& z,
                                                        const T& alpha) {
  const double z_dbl = value_of_rec(z);
  const double alpha_dbl = value_of_rec(alpha);
  if (z_dbl == NEGATIVE_INFTY) {
    return {T(NEGATIVE_INFTY), T(0.0), T(0.0)};
  }
  if (z_dbl == INFTY) {
    return {T(0.0), T(0.0), T(0.0)};
  }
  if (z_dbl < 0 && alpha_dbl > 0 && -alpha_dbl * z_dbl >= 2.5) {
    return skew_normal_log_cdf_tail(z, alpha);
  }
  const T alpha_z = alpha * z;
  const T cdf
      = alpha_dbl > 1
            ? T(2.0 * owens_t(abs(alpha_z), inv(alpha))
                + 0.5 * erfc(-alpha_z * INV_SQRT_TWO) * erf(z * INV_SQRT_TWO))
            : T(0.5 * erfc(-z * INV_SQRT_TWO) - 2.0 * owens_t(z, alpha));
  const T log_cdf = log(cdf);
  // exp(-z^2 / 2) and exp(-z^2 (1 + alpha^2) / 2) can be subnormal or 0
  // where their ratio to the cdf is not
  const T z_sq = square(z);
  const T d_z = exp(log(erfc(-alpha_z * INV_SQRT_TWO)) + NEG_LOG_SQRT_TWO_PI
                    - 0.5 * z_sq - log_cdf);
  const T d_alpha = -exp(-0.5 * (z_sq + square(alpha_z)) - LOG_PI
                         - log1p(square(alpha)) - log_cdf);
  return {log_cdf, d_z, d_alpha};
}

/**
 * Log of the skew normal cdf. The template argument `func` is the function
 * name in error messages, so that `skew_normal_lccdf` can call this
 * implementation with reflected arguments.
 *
 * @tparam func name of the calling function, for error messages
 * @tparam T_y type of the random variable
 * @tparam T_loc type of the location parameter
 * @tparam T_scale type of the scale parameter
 * @tparam T_shape type of the shape parameter
 * @param y random variable
 * @param mu location parameter
 * @param sigma scale parameter
 * @param alpha shape parameter
 * @return The log of the skew normal cdf evaluated at the specified
 *   arguments. If given containers, the log of the product of the cdfs.
 */
template <const char* func, typename T_y, typename T_loc, typename T_scale,
          typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lcdf_impl(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using T_alpha_ref = ref_type_t<T_shape>;
  static constexpr const char* function = func;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape parameter",
                         alpha);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_alpha_ref alpha_ref = alpha;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive(function, "Scale parameter", sigma_ref);
  check_finite(function, "Shape parameter", alpha_ref);

  if (size_zero(y, mu, sigma, alpha)) {
    return 0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, alpha_ref);
  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  const size_t N = max_size(y, mu, sigma, alpha);

  T_partials_return cdf_log(0.0);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return sigma_dbl = sigma_vec.val(n);
    const T_partials_return z = (y_vec.val(n) - mu_vec.val(n)) / sigma_dbl;
    const auto terms
        = skew_normal_log_cdf<T_partials_return>(z, alpha_vec.val(n));
    cdf_log += terms.log_cdf;
    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
      const T_partials_return d_y = terms.d_z / sigma_dbl;
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += d_y;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= d_y;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n] -= d_y * z;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<3>(ops_partials)[n] += terms.d_alpha;
    }
  }
  return ops_partials.build(cdf_log);
}
}  // namespace internal

/** \ingroup prob_dists
 * Returns the log of the skew normal cumulative distribution function.
 *
 * @tparam T_y type of the random variable
 * @tparam T_loc type of the location parameter
 * @tparam T_scale type of the scale parameter
 * @tparam T_shape type of the shape parameter
 * @param y random variable
 * @param mu location parameter
 * @param sigma scale parameter
 * @param alpha shape parameter
 * @return The log of the skew normal cdf evaluated at the specified
 *   arguments. If given containers, the log of the product of the cdfs.
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lcdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return internal::skew_normal_lcdf_impl<internal::skew_normal_lcdf_func>(
      y, mu, sigma, alpha);
}

}  // namespace math
}  // namespace stan
#endif
