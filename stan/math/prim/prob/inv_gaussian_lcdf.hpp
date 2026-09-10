#ifndef STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LCDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAUSSIAN_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/inv_square.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {
namespace internal {

// erfc underflows to zero below z = -38.5; the two branches agree to a
// relative 3e-15 at this cutoff.
static constexpr double LOG_PHI_ASYMPTOTIC_CUTOFF = -30.0;

/**
 * Return the log of the standard normal cumulative distribution function. At
 * or above LOG_PHI_ASYMPTOTIC_CUTOFF this is log(0.5 * erfc(-z / sqrt(2)));
 * below it the asymptotic expansion
 * \f$-z^2/2 - \log(-z) - \log(2\pi)/2
 * + \log(1 - s + 3s^2 - 15s^3 + 105s^4)\f$, \f$s = z^{-2}\f$, carries the
 * lower tail to a relative 1e-15 down to z = -1000.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T> log_Phi(const T& z) {
  using std::log;
  if (value_of_rec(z) >= LOG_PHI_ASYMPTOTIC_CUTOFF) {
    return LOG_HALF + log(erfc(-z * INV_SQRT_TWO));
  }
  const auto s = inv_square(z);
  return -0.5 * square(z) - log(-z) - HALF_LOG_TWO_PI
         + log1p(s * (-1.0 + s * (3.0 + s * (-15.0 + s * 105.0))));
}

struct log_Phi_fun {
  template <typename T>
  static inline auto fun(T&& z) {
    return log_Phi(std::forward<T>(z));
  }
};

/**
 * A vectorized version of log_Phi().
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto log_Phi(T&& z) {
  return apply_scalar_unary<log_Phi_fun, T>::apply(std::forward<T>(z));
}

/**
 * Return the log of the standard normal density at the specified value.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline return_type_t<T> log_std_normal_density(const T& z) {
  return -0.5 * square(z) - HALF_LOG_TWO_PI;
}

struct log_std_normal_density_fun {
  template <typename T>
  static inline auto fun(T&& z) {
    return log_std_normal_density(std::forward<T>(z));
  }
};

/**
 * A vectorized version of log_std_normal_density().
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto log_std_normal_density(T&& z) {
  return apply_scalar_unary<log_std_normal_density_fun, T>::apply(
      std::forward<T>(z));
}

/**
 * Return the log of the second inverse Gaussian CDF term,
 * \f$\log\left(e^{2\lambda/\mu}\,\Phi(-z_2)\right)\f$. Since
 * \f$z_2^2 - z_1^2 = 4\lambda/\mu\f$, the asymptotic expansion of
 * \f$\log\Phi(-z_2)\f$ cancels \f$2\lambda/\mu\f$ against \f$-z_2^2/2\f$ and
 * leaves \f$\log\phi(z_1) - \log z_2 + \log(1 - s + 3s^2 - 15s^3 + 105s^4)\f$
 * with \f$s = z_2^{-2}\f$, in which no two terms are large and opposed. Below
 * the asymptotic regime \f$z_2 \le 30\f$ bounds \f$2\lambda/\mu\f$ by 450 and
 * the direct form is accurate.
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_scaled_upper_term(const T1& z1, const T2& z2) {
  using std::log;
  if (value_of_rec(z2) > -LOG_PHI_ASYMPTOTIC_CUTOFF) {
    const auto s = inv_square(z2);
    return log_std_normal_density(z1) - log(z2)
           + log1p(s * (-1.0 + s * (3.0 + s * (-15.0 + s * 105.0))));
  }
  return 0.5 * (square(z2) - square(z1)) + log_Phi(-z2);
}

/**
 * A vectorized version of log_scaled_upper_term().
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_scaled_upper_term(T1&& z1, T2&& z2) {
  return apply_scalar_binary(
      [](auto&& a, auto&& b) {
        return log_scaled_upper_term(std::forward<decltype(a)>(a),
                                     std::forward<decltype(b)>(b));
      },
      std::forward<T1>(z1), std::forward<T2>(z2));
}

}  // namespace internal

/** \ingroup prob_dists
 * Returns the inverse Gaussian log cumulative distribution function for the
 * given random variable, mean and shape. Given containers of matching sizes,
 * returns the log of the product of probabilities.
 *
 * <p>The CDF is \f$F(y \mid \mu, \lambda) = \Phi(z_1)
 * + e^{2\lambda/\mu}\Phi(-z_2)\f$ with
 * \f$z_{1,2} = \sqrt{\lambda/y}\,(y/\mu \mp 1)\f$. The factor
 * \f$e^{2\lambda/\mu}\f$ overflows a double above \f$2\lambda/\mu = 710\f$,
 * so the two terms are combined in log space. That grouping of \f$z\f$ is
 * exact in floating point at \f$y = \mu\f$.
 *
 * <p>Both boundaries of the support are handled elementwise; the partials
 * are zero at both.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of mean parameter
 * @tparam T_shape type of shape parameter
 * @param y scalar or container of scalars
 * @param mu mean parameter
 * @param lambda shape parameter
 * @return the log of the product of probabilities
 * @throw std::domain_error if the mean or the shape is not positive and
 * finite, if the random variable is negative, or if any argument is NaN.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <typename T_y, typename T_loc, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_shape> inv_gaussian_lcdf(
    const T_y& y, const T_loc& mu, const T_shape& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "inv_gaussian_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Mean parameter", mu,
                         "Shape parameter", lambda);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Mean parameter", mu_val);
  check_positive_finite(function, "Shape parameter", lambda_val);

  if (size_zero(y, mu, lambda)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, lambda_ref);

  // Boundary masks. The standardized arguments are 0 * inf at y == inf, so
  // that contribution is substituted; at y == 0 they reach their limits on
  // their own.
  const auto& is_inf = to_ref(y_val == INFTY);
  const auto& is_bdry = to_ref((y_val == 0) || (y_val == INFTY));

  constexpr bool any_ad = is_any_autodiff_v<T_y, T_loc, T_shape>;

  const auto& inv_mu = to_ref_if<any_ad>(inv(mu_val));
  const auto& inv_y = to_ref_if<any_ad>(inv(y_val));
  const auto& sqrt_lambda_over_y = to_ref_if<any_ad>(sqrt(lambda_val * inv_y));
  const auto& y_over_mu = to_ref_if<any_ad>(y_val * inv_mu);
  const auto& z1 = to_ref(sqrt_lambda_over_y * (y_over_mu - 1.0));
  const auto& z2 = to_ref(sqrt_lambda_over_y * (y_over_mu + 1.0));

  const auto& log_upper = to_ref(internal::log_scaled_upper_term(z1, z2));
  const auto& lcdf_elt
      = to_ref(select(is_inf, T_partials_return(0),
                      log_sum_exp(internal::log_Phi(z1), log_upper)));

  T_partials_return lcdf = sum(lcdf_elt);

  if constexpr (any_ad) {
    // exp(2 lambda / mu) phi(z2) == phi(z1) exactly, because
    // z2^2 - z1^2 == 4 lambda / mu, which collapses both phi terms into one.
    //
    // Each partial is masked separately: at y == inf the shape partial is
    // 0 / 0 on its own. The inner select covers elements whose log CDF has
    // saturated to -inf at an interior y.
    const auto& is_underflow = to_ref(lcdf_elt == NEGATIVE_INFTY);
    const auto& w_dens
        = to_ref(exp(internal::log_std_normal_density(z1) - lcdf_elt));
    const auto& w_upper = to_ref(exp(log_upper - lcdf_elt));
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials)
          = select(is_bdry, T_partials_return(0),
                   select(is_underflow, T_partials_return(0),
                          w_dens * sqrt_lambda_over_y * inv_y));
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials)
          = select(is_bdry, T_partials_return(0),
                   select(is_underflow, T_partials_return(0),
                          -2.0 * lambda_val * w_upper * square(inv_mu)));
    }
    if constexpr (is_autodiff_v<T_shape>) {
      partials<2>(ops_partials) = select(
          is_bdry, T_partials_return(0),
          select(is_underflow, T_partials_return(0),
                 2.0 * w_upper * inv_mu - w_dens * inv_y / sqrt_lambda_over_y));
    }
  }
  return ops_partials.build(lcdf);
}

}  // namespace math
}  // namespace stan
#endif
