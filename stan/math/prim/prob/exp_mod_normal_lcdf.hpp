#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/exp_mod_normal_check_args.hpp>
#include <array>
#include <cstddef>
#include <utility>

namespace stan {
namespace math {
namespace internal {

constexpr double exp_mod_normal_large_p = 0.75;
constexpr double exp_mod_normal_log_large_p = -0.2876820724517809;
constexpr double exp_mod_normal_b_only_t = -6.5;
constexpr double exp_mod_normal_erfcx_close_ratio = 0.99;
constexpr double exp_mod_normal_erfcx_tail_start = 4.0;

/**
 * Returns `h` such that
 * `erfcx(lo) - erfcx(hi) = erfcx(lo) * d * (1 - lo * h) / hi` for
 * `hi = lo + d`, from the rational tail approximation of `erfcx`. Valid for
 * `lo >= exp_mod_normal_erfcx_tail_start`.
 *
 * @tparam T1 type of `lo`
 * @tparam T2 type of `hi`
 * @param lo lower argument
 * @param hi upper argument
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> exp_mod_normal_erfcx_tail_h_per_d(const T1& lo,
                                                               const T2& hi) {
  using T_return = return_type_t<T1, T2>;
  const auto w0 = inv(lo * lo);
  const auto w1 = inv(hi * hi);
  T_return a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
  T_return pw0 = 1, pw1 = 1, sum_k = 0;
  for (std::size_t k = 0; k < erfcx_tail_p.size(); ++k) {
    q0 += erfcx_tail_q[k] * pw0;
    q1 += erfcx_tail_q[k] * pw1;
    dd_q += erfcx_tail_q[k] * sum_k;
    sum_k = sum_k * w0 + pw1;
    a0 += erfcx_tail_p[k] * pw0 * w0;
    dd_a += erfcx_tail_p[k] * sum_k;
    pw0 *= w0;
    pw1 *= w1;
  }
  const auto pq = erfcx_tail_polynomials<true>(w0);
  return (SQRT_PI * (-w0 * ((hi + lo) * w1)) * (dd_a * q0 - a0 * dd_q)
          / (q1 * q0) / (1 + SQRT_PI * w0 * pq.first / pq.second));
}

/**
 * Returns `(erfcx(mid - d / 2) - erfcx(mid + d / 2)) / d` from the Taylor
 * series of `erfcx` around `mid`.
 *
 * @tparam T1 type of `mid`
 * @tparam T2 type of `d`
 * @param mid midpoint of the two arguments
 * @param d distance between the two arguments
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> exp_mod_normal_erfcx_mid_sum(const T1& mid,
                                                          const T2& d) {
  // e[k] is the k-th derivative of erfcx at mid
  std::array<return_type_t<T1, T2>, 8> e;
  e[0] = erfcx(mid);
  e[1] = erfcx_derivative(mid, e[0]);
  for (std::size_t k = 1; k < 7; ++k) {
    e[k + 1] = 2 * mid * e[k] + 2 * k * e[k - 1];
  }
  const auto h = 0.5 * d;
  const auto h2 = h * h;
  return -(e[1] + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
}

template <bool upper, typename T_y, typename T_loc, typename T_scale,
          typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lcdf_impl(
    const char* function, T_y&& y, T_loc&& mu, T_scale&& sigma,
    T_inv_scale&& lambda) {
  using T_partials = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  T_lambda_ref lambda_ref = std::forward<T_inv_scale>(lambda);
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));
  const bool y_inf = exp_mod_normal_check_args(function, y_val, mu_val,
                                               sigma_val, lambda_val);

  if (size_zero(y_ref, mu_ref, sigma_ref, lambda_ref)) {
    // empty input: the sum has no terms
    return 0.0;
  } else if (unlikely(y_inf)) {
    // an infinite y decides the result without evaluating any term
    if (any(y_val == NEGATIVE_INFTY)) {
      // y = -inf: the cdf is 0 and the ccdf is 1
      return upper ? 0.0 : NEGATIVE_INFTY;
    }
    return upper ? NEGATIVE_INFTY : 0.0;
  }

  scalar_seq_view<std::decay_t<decltype(y_val)>> y_vec(y_val);
  scalar_seq_view<std::decay_t<decltype(mu_val)>> mu_vec(mu_val);
  scalar_seq_view<std::decay_t<decltype(sigma_val)>> sigma_vec(sigma_val);
  scalar_seq_view<std::decay_t<decltype(lambda_val)>> lambda_vec(lambda_val);
  const std::size_t N = max_size(y_ref, mu_ref, sigma_ref, lambda_ref);
  T_partials lp = 0.0;
  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);
  for (std::size_t n = 0; n < N; ++n) {
    const auto sigma_n = sigma_vec[n];
    const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
    const auto v = lambda_vec[n] * sigma_n;
    const auto u = -z * INV_SQRT_TWO;
    const auto d = v * INV_SQRT_TWO;
    const auto t = u + d;
    const auto log_g = LOG_HALF - (0.5 * z) * z;
    if constexpr (upper) {
      if (z <= 0) {
        // y at or below mu: u and t are nonnegative, so erfcx cannot overflow
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        const auto g = exp(log_g);
        T_partials f;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          // erfcx(u) and erfcx(t) differ enough to subtract directly
          f = g * (e_u - e_t);
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          // close values in the erfcx tail: difference of the rational form
          f = g
              * (e_u
                 * (d * (1 - u * exp_mod_normal_erfcx_tail_h_per_d(u, t)) / t));
        } else {
          // close values below the tail: Taylor series around the midpoint
          f = g * (d * exp_mod_normal_erfcx_mid_sum(u + 0.5 * d, d));
        }
        const auto scale = -g / (1 - f);
        lp += log1p(-f);
        if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
          const auto d_z = (scale * v * e_t) / sigma_n;
          if constexpr (is_autodiff_v<T_y>) {
            partials<0>(ops_partials)[n] += d_z;
          }
          if constexpr (is_autodiff_v<T_loc>) {
            partials<1>(ops_partials)[n] -= d_z;
          }
        }
        if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
          const auto de_t = erfcx_derivative(t, e_t);
          if constexpr (is_autodiff_v<T_scale>) {
            partials<2>(ops_partials)[n]
                += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t)) / sigma_n;
          }
          if constexpr (is_autodiff_v<T_inv_scale>) {
            partials<3>(ops_partials)[n]
                += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
          }
        }
      } else if (const auto log_b = d * (2 * u + d);
                 t < exp_mod_normal_b_only_t
                 && log_b < exp_mod_normal_log_large_p) {
        // far right tail: the erfcx terms are below rounding, ccdf = exp(log_b)
        lp += log_b;
        if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
          const auto d_z = -v / sigma_n;
          if constexpr (is_autodiff_v<T_y>) {
            partials<0>(ops_partials)[n] += d_z;
          }
          if constexpr (is_autodiff_v<T_loc>) {
            partials<1>(ops_partials)[n] -= d_z;
          }
        }
        if constexpr (is_autodiff_v<T_scale>) {
          partials<2>(ops_partials)[n] += (v * v) / sigma_n;
        }
        if constexpr (is_autodiff_v<T_inv_scale>) {
          partials<3>(ops_partials)[n] += (v - z) * sigma_n;
        }
      } else if (t >= 0) {
        // y above mu, t >= 0: erfcx(-u) and erfcx(t) cannot overflow
        const auto g = exp(log_g);
        const auto e_t = erfcx(t);
        const auto e_sum = erfcx(-u) + e_t;
        if (g * e_sum > exp_mod_normal_large_p) {
          // ccdf near 1: its log is log1p(-cdf), with the cdf computed directly
          T_partials f;
          if (const auto e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            // erfcx(u) and erfcx(t) differ enough to subtract directly
            f = g * (e_u - e_t);
          } else {
            // close values (u < 0 here): Taylor series around the midpoint
            f = g * (d * exp_mod_normal_erfcx_mid_sum(u + 0.5 * d, d));
          }
          const auto scale = -g / (1 - f);
          lp += log1p(-f);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (scale * v * e_t) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto de_t = erfcx_derivative(t, e_t);
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t)) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
            }
          }
        } else {
          // ccdf away from 1: take its log directly
          lp += log_g + log(e_sum);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (-v * e_t / e_sum) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto de_t = erfcx_derivative(t, e_t);
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (-v * SQRT_TWO * (u * e_t - 0.5 * de_t) / e_sum) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += (de_t * INV_SQRT_TWO / e_sum) * sigma_n;
            }
          }
        }
      } else {
        // t < 0: erfcx(t) can overflow, so use erfcx(-t) and exp(log_b)
        const auto e_mu = erfcx(-u);
        const auto e_mt = erfcx(-t);
        const auto half_exp_mt2 = 0.5 * exp(-t * t);
        const auto rel_tail = half_exp_mt2 * (e_mt - e_mu);
        const auto p_over_b = 1 - rel_tail;
        const auto exp_log_b = exp(log_b);
        if (exp_log_b * p_over_b > exp_mod_normal_large_p) {
          // ccdf near 1: its log is log1p(-cdf), with the cdf computed directly
          const auto g = exp_log_b * half_exp_mt2;
          const auto b = exp_log_b * (1 - half_exp_mt2 * e_mt);
          T_partials f;
          if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
            // erfcx(-t) and erfcx(-u) differ enough to subtract directly
            f = -expm1(log_b) + g * (e_mt - e_mu);
          } else if (const auto lo = -t;
                     lo >= exp_mod_normal_erfcx_tail_start) {
            // close values in the erfcx tail: difference of the rational form
            const auto hi = -u;
            f = -expm1(log_b)
                + g
                      * (e_mt
                         * (d
                            * (1
                               - lo * exp_mod_normal_erfcx_tail_h_per_d(lo, hi))
                            / hi));
          } else {
            // close values below the tail: Taylor series around the midpoint
            f = -expm1(log_b)
                + g * (d * exp_mod_normal_erfcx_mid_sum(-u - 0.5 * d, d));
          }
          const auto divisor = f - 1;
          lp += log1p(-f);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (v * b / divisor) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (v * (phi - v * b) / divisor) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += ((phi - (v - z) * b) / divisor) * sigma_n;
            }
          }
        } else {
          // ccdf away from 1: log of exp(log_b) times its relative correction
          const auto b_p = (1 - half_exp_mt2 * e_mt) / p_over_b;
          lp += log_b + log1p(-rel_tail);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (-v * b_p) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto phi_p = half_exp_mt2 * SQRT_TWO_OVER_SQRT_PI / p_over_b;
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (-v * (phi_p - v * b_p)) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += (-(phi_p - (v - z) * b_p)) * sigma_n;
            }
          }
        }
      }
    } else {
      if (z <= 0) {
        // y at or below mu: u and t are nonnegative, so erfcx cannot overflow
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        T_partials log_delta;
        T_partials d_e_t_over_delta;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          // erfcx(u) and erfcx(t) differ enough to subtract directly
          const auto delta = e_u - e_t;
          log_delta = log(delta);
          d_e_t_over_delta = d * (e_t / delta);
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          // close values in the erfcx tail: difference of the rational form
          const auto h_per_d = exp_mod_normal_erfcx_tail_h_per_d(u, t);
          const auto one_minus_u_h = 1 - u * h_per_d;
          log_delta = log(e_u) + log(d) + log(one_minus_u_h / t);
          d_e_t_over_delta = u * (1 + d * h_per_d) / one_minus_u_h;
        } else {
          // close values below the tail: Taylor series around the midpoint
          const auto neg_s = exp_mod_normal_erfcx_mid_sum(u + 0.5 * d, d);
          log_delta = log(d) + log(neg_s);
          d_e_t_over_delta = e_t / neg_s;
        }
        lp += log_g + log_delta;
        if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
          const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
          if constexpr (is_autodiff_v<T_y>) {
            partials<0>(ops_partials)[n] += d_z;
          }
          if constexpr (is_autodiff_v<T_loc>) {
            partials<1>(ops_partials)[n] -= d_z;
          }
        }
        if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
          const auto s = erfcx_derivative(t, e_t) / e_t;
          if constexpr (is_autodiff_v<T_scale>) {
            partials<2>(ops_partials)[n]
                += (d_e_t_over_delta * (2 * u - s)) / sigma_n;
          }
          if constexpr (is_autodiff_v<T_inv_scale>) {
            partials<3>(ops_partials)[n]
                += (-(s * d_e_t_over_delta) / v) * sigma_n;
          }
        }
      } else if (const auto log_b = d * (2 * u + d);
                 t < exp_mod_normal_b_only_t
                 && log_b < exp_mod_normal_log_large_p) {
        // far right tail: the erfcx terms are below rounding, ccdf = exp(log_b)
        const auto b = exp(log_b);
        const auto r = b / (1 - b);
        lp += log1p(-b);
        if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
          const auto d_z = (v * r) / sigma_n;
          if constexpr (is_autodiff_v<T_y>) {
            partials<0>(ops_partials)[n] += d_z;
          }
          if constexpr (is_autodiff_v<T_loc>) {
            partials<1>(ops_partials)[n] -= d_z;
          }
        }
        if constexpr (is_autodiff_v<T_scale>) {
          partials<2>(ops_partials)[n] += (-v * v * r) / sigma_n;
        }
        if constexpr (is_autodiff_v<T_inv_scale>) {
          partials<3>(ops_partials)[n] += ((z - v) * r) * sigma_n;
        }
      } else if (t >= 0) {
        // y above mu, t >= 0: erfcx(-u) and erfcx(t) cannot overflow
        const auto g = exp(log_g);
        const auto e_t = erfcx(t);
        const auto b = g * e_t;
        if (const auto p = g * erfcx(-u) + b; p > exp_mod_normal_large_p) {
          // cdf small: 1 - p would cancel, so compute the cdf directly
          T_partials log_delta;
          T_partials d_e_t_over_delta;
          if (const auto e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            // erfcx(u) and erfcx(t) differ enough to subtract directly
            const auto delta = e_u - e_t;
            log_delta = log(delta);
            d_e_t_over_delta = d * (e_t / delta);
          } else {
            // close values (u < 0 here): Taylor series around the midpoint
            const auto neg_s = exp_mod_normal_erfcx_mid_sum(u + 0.5 * d, d);
            log_delta = log(d) + log(neg_s);
            d_e_t_over_delta = e_t / neg_s;
          }
          lp += log_g + log_delta;
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto s = erfcx_derivative(t, e_t) / e_t;
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (d_e_t_over_delta * (2 * u - s)) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += (-(s * d_e_t_over_delta) / v) * sigma_n;
            }
          }
        } else {
          // cdf away from 0: 1 - p is accurate
          const auto f = 1 - p;
          lp += log1p(-p);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (v * b / f) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto de_t = erfcx_derivative(t, e_t);
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n]
                  += (v * g * SQRT_TWO * (u * e_t - 0.5 * de_t) / f) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += (-g * de_t * INV_SQRT_TWO / f) * sigma_n;
            }
          }
        }
      } else {
        // t < 0: erfcx(t) can overflow, so use erfcx(-t) and exp(log_b)
        const auto g = exp(log_g);
        const auto e_mu = erfcx(-u);
        const auto e_mt = erfcx(-t);
        const auto b = exp(log_b) - g * e_mt;
        if (const auto p = g * e_mu + b; p > exp_mod_normal_large_p) {
          // cdf small: 1 - p would cancel, so compute the cdf directly
          T_partials f;
          if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
            // erfcx(-t) and erfcx(-u) differ enough to subtract directly
            f = -expm1(log_b) + g * (e_mt - e_mu);
          } else if (const auto lo = -t;
                     lo >= exp_mod_normal_erfcx_tail_start) {
            // close values in the erfcx tail: difference of the rational form
            const auto hi = -u;
            f = -expm1(log_b)
                + g
                      * (e_mt
                         * (d
                            * (1
                               - lo * exp_mod_normal_erfcx_tail_h_per_d(lo, hi))
                            / hi));
          } else {
            // close values below the tail: Taylor series around the midpoint
            f = -expm1(log_b)
                + g * (d * exp_mod_normal_erfcx_mid_sum(-u - 0.5 * d, d));
          }
          lp += log(f);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (v * b / f) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n] += (v * (phi - v * b) / f) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += ((phi - (v - z) * b) / f) * sigma_n;
            }
          }
        } else {
          // cdf away from 0: 1 - p is accurate
          const auto f = 1 - p;
          lp += log1p(-p);
          if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
            const auto d_z = (v * b / f) / sigma_n;
            if constexpr (is_autodiff_v<T_y>) {
              partials<0>(ops_partials)[n] += d_z;
            }
            if constexpr (is_autodiff_v<T_loc>) {
              partials<1>(ops_partials)[n] -= d_z;
            }
          }
          if constexpr (is_autodiff_v<T_scale> || is_autodiff_v<T_inv_scale>) {
            const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
            if constexpr (is_autodiff_v<T_scale>) {
              partials<2>(ops_partials)[n] += (v * (phi - v * b) / f) / sigma_n;
            }
            if constexpr (is_autodiff_v<T_inv_scale>) {
              partials<3>(ops_partials)[n]
                  += ((phi - (v - z) * b) / f) * sigma_n;
            }
          }
        }
      }
    }
  }
  return ops_partials.build(lp);
}

}  // namespace internal

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lcdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  return internal::exp_mod_normal_lcdf_impl<false>(
      "exp_mod_normal_lcdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma), std::forward<T_inv_scale>(lambda));
}

}  // namespace math
}  // namespace stan
#endif
