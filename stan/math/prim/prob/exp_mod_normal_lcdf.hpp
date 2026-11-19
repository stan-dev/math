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
#include <stan/math/prim/fun/log_sum_exp.hpp>
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
  constexpr bool any_autodiff
      = is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>;
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
    return 0.0;
  } else if (unlikely(y_inf)) {
    if (any(y_val == NEGATIVE_INFTY)) {
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
  if constexpr (upper) {
    if constexpr (any_autodiff) {
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
        if (z <= 0) {
          if (const auto e_t = erfcx(t), e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            const auto g = exp(log_g);
            const auto f = g * (e_u - e_t);
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto de_t = erfcx_derivative(t, e_t);
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t))
                       / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
              }
            }
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
            T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
            const auto g = exp(log_g);
            const auto f
                = g
                  * (e_u
                     * (d
                        * (1
                           - u
                                 * (SQRT_PI * (-w0 * ((t + u) * w1))
                                    * (dd_a * q0 - a0 * dd_q) / (q1 * q0)
                                    / (1
                                       + SQRT_PI * w0 * pq.first / pq.second)))
                        / t));
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto de_t = erfcx_derivative(t, e_t);
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t))
                       / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
              }
            }
          } else {
            const auto mid = u + 0.5 * d;
            T_partials e_km1 = erfcx(mid);
            T_partials e_k = erfcx_derivative(mid, e_km1);
            std::array<T_partials, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto g = exp(log_g);
            const auto f
                = g
                  * (d
                     * -(odd[0]
                         + h2
                               * (odd[1] / 6
                                  + h2 * (odd[2] / 120 + h2 * odd[3] / 5040))));
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto de_t = erfcx_derivative(t, e_t);
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t))
                       / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
              }
            }
          }
        } else if (const auto log_b = d * (2 * u + d);
                   t < exp_mod_normal_b_only_t
                   && log_b < exp_mod_normal_log_large_p) {
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
          const auto g = exp(log_g);
          const auto e_t = erfcx(t);
          const auto e_sum = erfcx(-u) + e_t;
          if (g * e_sum > exp_mod_normal_large_p) {
            if (const auto e_u = erfcx(u);
                e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
              const auto f = g * (e_u - e_t);
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                const auto de_t = erfcx_derivative(t, e_t);
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t))
                         / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
                }
              }
            } else {
              const auto mid = u + 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              const auto f = g
                             * (d
                                * -(odd[0]
                                    + h2
                                          * (odd[1] / 6
                                             + h2
                                                   * (odd[2] / 120
                                                      + h2 * odd[3] / 5040))));
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                const auto de_t = erfcx_derivative(t, e_t);
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t))
                         / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += (-scale * de_t * INV_SQRT_TWO) * sigma_n;
                }
              }
            }
          } else {
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto de_t = erfcx_derivative(t, e_t);
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (-v * SQRT_TWO * (u * e_t - 0.5 * de_t) / e_sum)
                       / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (de_t * INV_SQRT_TWO / e_sum) * sigma_n;
              }
            }
          }
        } else {
          const auto e_mu = erfcx(-u);
          const auto e_mt = erfcx(-t);
          const auto half_exp_mt2 = 0.5 * exp(-t * t);
          const auto rel_tail = half_exp_mt2 * (e_mt - e_mu);
          const auto p_over_b = 1 - rel_tail;
          const auto exp_log_b = exp(log_b);
          if (exp_log_b * p_over_b > exp_mod_normal_large_p) {
            if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
              const auto g = exp_log_b * half_exp_mt2;
              const auto b = exp_log_b * (1 - half_exp_mt2 * e_mt);
              const auto f = -expm1(log_b) + g * (e_mt - e_mu);
              const auto divisor = -(1 - f);
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
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
            } else if (const auto lo = -t;
                       lo >= exp_mod_normal_erfcx_tail_start) {
              const auto hi = -u;
              const auto w0 = inv(lo * lo);
              const auto w1 = inv(hi * hi);
              T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
              T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
              const auto g = exp_log_b * half_exp_mt2;
              const auto b = exp_log_b * (1 - half_exp_mt2 * e_mt);
              const auto f
                  = -expm1(log_b)
                    + g
                          * (e_mt
                             * (d
                                * (1
                                   - lo
                                         * (SQRT_PI * (-w0 * ((hi + lo) * w1))
                                            * (dd_a * q0 - a0 * dd_q)
                                            / (q1 * q0)
                                            / (1
                                               + SQRT_PI * w0 * pq.first
                                                     / pq.second)))
                                / hi));
              const auto divisor = -(1 - f);
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
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
              const auto mid = -u - 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              const auto g = exp_log_b * half_exp_mt2;
              const auto b = exp_log_b * (1 - half_exp_mt2 * e_mt);
              const auto f
                  = -expm1(log_b)
                    + g
                          * (d
                             * -(odd[0]
                                 + h2
                                       * (odd[1] / 6
                                          + h2
                                                * (odd[2] / 120
                                                   + h2 * odd[3] / 5040))));
              const auto divisor = -(1 - f);
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
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
            }
          } else {
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto phi_p
                  = half_exp_mt2 * SQRT_TWO_OVER_SQRT_PI / p_over_b;
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
      }
      return ops_partials.build(lp);
    } else {
      for (std::size_t n = 0; n < N; ++n) {
        const auto sigma_n = sigma_vec[n];
        const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
        const auto u = -z * INV_SQRT_TWO;
        const auto d = (lambda_vec[n] * sigma_n) * INV_SQRT_TWO;
        const auto t = u + d;
        const auto log_g = LOG_HALF - (0.5 * z) * z;
        if (z <= 0) {
          if (const auto e_t = erfcx(t), e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            lp += log1p(-(exp(log_g) * (e_u - e_t)));
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
            T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
            lp += log1p(-(
                exp(log_g)
                * (e_u
                   * (d
                      * (1
                         - u
                               * (SQRT_PI * (-w0 * ((t + u) * w1))
                                  * (dd_a * q0 - a0 * dd_q) / (q1 * q0)
                                  / (1 + SQRT_PI * w0 * pq.first / pq.second)))
                      / t))));
          } else {
            const auto mid = u + 0.5 * d;
            T_partials e_km1 = erfcx(mid);
            T_partials e_k = erfcx_derivative(mid, e_km1);
            std::array<T_partials, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            lp += log1p(-(
                exp(log_g)
                * (d
                   * -(odd[0]
                       + h2
                             * (odd[1] / 6
                                + h2 * (odd[2] / 120 + h2 * odd[3] / 5040))))));
          }
        } else if (const auto log_b = d * (2 * u + d);
                   t < exp_mod_normal_b_only_t
                   && log_b < exp_mod_normal_log_large_p) {
          lp += log_b;
        } else if (t >= 0) {
          const auto g = exp(log_g);
          const auto e_t = erfcx(t);
          const auto e_sum = erfcx(-u) + e_t;
          if (g * e_sum > exp_mod_normal_large_p) {
            if (const auto e_u = erfcx(u);
                e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
              lp += log1p(-(g * (e_u - e_t)));
            } else {
              const auto mid = u + 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              lp += log1p(-(g
                            * (d
                               * -(odd[0]
                                   + h2
                                         * (odd[1] / 6
                                            + h2
                                                  * (odd[2] / 120
                                                     + h2 * odd[3] / 5040))))));
            }
          } else {
            lp += log_g + log(e_sum);
          }
        } else {
          const auto e_mu = erfcx(-u);
          const auto e_mt = erfcx(-t);
          const auto half_exp_mt2 = 0.5 * exp(-t * t);
          const auto rel_tail = half_exp_mt2 * (e_mt - e_mu);
          const auto exp_log_b = exp(log_b);
          if (exp_log_b * (1 - rel_tail) > exp_mod_normal_large_p) {
            if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
              lp += log1p(-(-expm1(log_b)
                            + (exp_log_b * half_exp_mt2) * (e_mt - e_mu)));
            } else if (const auto lo = -t;
                       lo >= exp_mod_normal_erfcx_tail_start) {
              const auto hi = -u;
              const auto w0 = inv(lo * lo);
              const auto w1 = inv(hi * hi);
              T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
              T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
              lp += log1p(-(
                  -expm1(log_b)
                  + (exp_log_b * half_exp_mt2)
                        * (e_mt
                           * (d
                              * (1
                                 - lo
                                       * (SQRT_PI * (-w0 * ((hi + lo) * w1))
                                          * (dd_a * q0 - a0 * dd_q) / (q1 * q0)
                                          / (1
                                             + SQRT_PI * w0 * pq.first
                                                   / pq.second)))
                              / hi))));
            } else {
              const auto mid = -u - 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              lp += log1p(
                  -(-expm1(log_b)
                    + (exp_log_b * half_exp_mt2)
                          * (d
                             * -(odd[0]
                                 + h2
                                       * (odd[1] / 6
                                          + h2
                                                * (odd[2] / 120
                                                   + h2 * odd[3] / 5040))))));
            }
          } else {
            lp += log_b + log1p(-rel_tail);
          }
        }
      }
      return lp;
    }
  } else {
    if constexpr (any_autodiff) {
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
        if (z <= 0) {
          if (const auto e_t = erfcx(t), e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            const auto delta = e_u - e_t;
            const auto d_e_t_over_delta = d * (e_t / delta);
            lp += log_g + log(delta);
            if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
              const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
              if constexpr (is_autodiff_v<T_y>) {
                partials<0>(ops_partials)[n] += d_z;
              }
              if constexpr (is_autodiff_v<T_loc>) {
                partials<1>(ops_partials)[n] -= d_z;
              }
            }
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              T_partials s;
              if (t > 1e150) {
                s = -inv(t);
              } else {
                s = erfcx_derivative(t, e_t) / e_t;
              }
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (d_e_t_over_delta * (2 * u - s)) / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-(s * d_e_t_over_delta) / v) * sigma_n;
              }
            }
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
            T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
            const auto h_per_d
                = (SQRT_PI * (-w0 * ((t + u) * w1)) * (dd_a * q0 - a0 * dd_q)
                   / (q1 * q0) / (1 + SQRT_PI * w0 * pq.first / pq.second));
            const auto one_minus_u_h = 1 - u * h_per_d;
            const auto d_e_t_over_delta = u * (1 + d * h_per_d) / one_minus_u_h;
            lp += log_g + (log(e_u) + log(d) + log(one_minus_u_h / t));
            if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
              const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
              if constexpr (is_autodiff_v<T_y>) {
                partials<0>(ops_partials)[n] += d_z;
              }
              if constexpr (is_autodiff_v<T_loc>) {
                partials<1>(ops_partials)[n] -= d_z;
              }
            }
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              T_partials s;
              if (t > 1e150) {
                s = -inv(t);
              } else {
                s = erfcx_derivative(t, e_t) / e_t;
              }
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
            const auto mid = u + 0.5 * d;
            T_partials e_km1 = erfcx(mid);
            T_partials e_k = erfcx_derivative(mid, e_km1);
            std::array<T_partials, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                odd[0]
                + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
            const auto d_e_t_over_delta = e_t / neg_s;
            lp += log_g + (log(d) + log(neg_s));
            if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
              const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
              if constexpr (is_autodiff_v<T_y>) {
                partials<0>(ops_partials)[n] += d_z;
              }
              if constexpr (is_autodiff_v<T_loc>) {
                partials<1>(ops_partials)[n] -= d_z;
              }
            }
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              T_partials s;
              if (t > 1e150) {
                s = -inv(t);
              } else {
                s = erfcx_derivative(t, e_t) / e_t;
              }
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (d_e_t_over_delta * (2 * u - s)) / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-(s * d_e_t_over_delta) / v) * sigma_n;
              }
            }
          }
        } else if (const auto log_b = d * (2 * u + d);
                   t < exp_mod_normal_b_only_t
                   && log_b < exp_mod_normal_log_large_p) {
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
          const auto g = exp(log_g);
          const auto e_t = erfcx(t);
          const auto b = g * e_t;
          if (const auto p = g * erfcx(-u) + b; p > exp_mod_normal_large_p) {
            if (const auto e_u = erfcx(u);
                e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
              const auto delta = e_u - e_t;
              const auto d_e_t_over_delta = d * (e_t / delta);
              lp += log_g + log(delta);
              if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
                const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
                if constexpr (is_autodiff_v<T_y>) {
                  partials<0>(ops_partials)[n] += d_z;
                }
                if constexpr (is_autodiff_v<T_loc>) {
                  partials<1>(ops_partials)[n] -= d_z;
                }
              }
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                T_partials s;
                if (t > 1e150) {
                  s = -inv(t);
                } else {
                  s = erfcx_derivative(t, e_t) / e_t;
                }
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
              const auto mid = u + 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              const auto neg_s
                  = -(odd[0]
                      + h2
                            * (odd[1] / 6
                               + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
              const auto d_e_t_over_delta = e_t / neg_s;
              lp += log_g + (log(d) + log(neg_s));
              if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
                const auto d_z = (SQRT_TWO * d_e_t_over_delta) / sigma_n;
                if constexpr (is_autodiff_v<T_y>) {
                  partials<0>(ops_partials)[n] += d_z;
                }
                if constexpr (is_autodiff_v<T_loc>) {
                  partials<1>(ops_partials)[n] -= d_z;
                }
              }
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                T_partials s;
                if (t > 1e150) {
                  s = -inv(t);
                } else {
                  s = erfcx_derivative(t, e_t) / e_t;
                }
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (d_e_t_over_delta * (2 * u - s)) / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += (-(s * d_e_t_over_delta) / v) * sigma_n;
                }
              }
            }
          } else {
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto de_t = erfcx_derivative(t, e_t);
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (v * g * SQRT_TWO * (u * e_t - 0.5 * de_t) / f)
                       / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += (-g * de_t * INV_SQRT_TWO / f) * sigma_n;
              }
            }
          }
        } else {
          const auto g = exp(log_g);
          const auto e_mu = erfcx(-u);
          const auto e_mt = erfcx(-t);
          const auto b = exp(log_b) - g * e_mt;
          if (const auto p = g * e_mu + b; p > exp_mod_normal_large_p) {
            if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
              const auto f = -expm1(log_b) + g * (e_mt - e_mu);
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (v * (phi - v * b) / f) / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += ((phi - (v - z) * b) / f) * sigma_n;
                }
              }
            } else if (const auto lo = -t;
                       lo >= exp_mod_normal_erfcx_tail_start) {
              const auto hi = -u;
              const auto w0 = inv(lo * lo);
              const auto w1 = inv(hi * hi);
              T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
              T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
              const auto f
                  = -expm1(log_b)
                    + g
                          * (e_mt
                             * (d
                                * (1
                                   - lo
                                         * (SQRT_PI * (-w0 * ((hi + lo) * w1))
                                            * (dd_a * q0 - a0 * dd_q)
                                            / (q1 * q0)
                                            / (1
                                               + SQRT_PI * w0 * pq.first
                                                     / pq.second)))
                                / hi));
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (v * (phi - v * b) / f) / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += ((phi - (v - z) * b) / f) * sigma_n;
                }
              }
            } else {
              const auto mid = -u - 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              const auto f
                  = -expm1(log_b)
                    + g
                          * (d
                             * -(odd[0]
                                 + h2
                                       * (odd[1] / 6
                                          + h2
                                                * (odd[2] / 120
                                                   + h2 * odd[3] / 5040))));
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
              if constexpr (is_autodiff_v<
                                T_scale> || is_autodiff_v<T_inv_scale>) {
                const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
                if constexpr (is_autodiff_v<T_scale>) {
                  partials<2>(ops_partials)[n]
                      += (v * (phi - v * b) / f) / sigma_n;
                }
                if constexpr (is_autodiff_v<T_inv_scale>) {
                  partials<3>(ops_partials)[n]
                      += ((phi - (v - z) * b) / f) * sigma_n;
                }
              }
            }
          } else {
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
            if constexpr (is_autodiff_v<
                              T_scale> || is_autodiff_v<T_inv_scale>) {
              const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
              if constexpr (is_autodiff_v<T_scale>) {
                partials<2>(ops_partials)[n]
                    += (v * (phi - v * b) / f) / sigma_n;
              }
              if constexpr (is_autodiff_v<T_inv_scale>) {
                partials<3>(ops_partials)[n]
                    += ((phi - (v - z) * b) / f) * sigma_n;
              }
            }
          }
        }
      }
      return ops_partials.build(lp);
    } else {
      for (std::size_t n = 0; n < N; ++n) {
        const auto sigma_n = sigma_vec[n];
        const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
        const auto u = -z * INV_SQRT_TWO;
        const auto d = (lambda_vec[n] * sigma_n) * INV_SQRT_TWO;
        const auto t = u + d;
        const auto log_g = LOG_HALF - (0.5 * z) * z;
        if (z <= 0) {
          if (const auto e_t = erfcx(t), e_u = erfcx(u);
              e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            lp += log_g + log(e_u - e_t);
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
            T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
            lp += log_g
                  + (log(e_u) + log(d)
                     + log(
                         (1
                          - u
                                * (SQRT_PI * (-w0 * ((t + u) * w1))
                                   * (dd_a * q0 - a0 * dd_q) / (q1 * q0)
                                   / (1 + SQRT_PI * w0 * pq.first / pq.second)))
                         / t));
          } else {
            const auto mid = u + 0.5 * d;
            T_partials e_km1 = erfcx(mid);
            T_partials e_k = erfcx_derivative(mid, e_km1);
            std::array<T_partials, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            lp += log_g
                  + (log(d)
                     + log(-(odd[0]
                             + h2
                                   * (odd[1] / 6
                                      + h2
                                            * (odd[2] / 120
                                               + h2 * odd[3] / 5040)))));
          }
        } else if (const auto log_b = d * (2 * u + d);
                   t < exp_mod_normal_b_only_t
                   && log_b < exp_mod_normal_log_large_p) {
          lp += log1p(-exp(log_b));
        } else if (t >= 0) {
          const auto g = exp(log_g);
          const auto e_t = erfcx(t);
          if (const auto p = g * erfcx(-u) + g * e_t;
              p > exp_mod_normal_large_p) {
            if (const auto e_u = erfcx(u);
                e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
              lp += log_g + log(e_u - e_t);
            } else {
              const auto mid = u + 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              lp += log_g
                    + (log(d)
                       + log(-(odd[0]
                               + h2
                                     * (odd[1] / 6
                                        + h2
                                              * (odd[2] / 120
                                                 + h2 * odd[3] / 5040)))));
            }
          } else {
            lp += log1p(-p);
          }
        } else {
          const auto g = exp(log_g);
          const auto e_mu = erfcx(-u);
          const auto e_mt = erfcx(-t);
          if (const auto p = g * e_mu + (exp(log_b) - g * e_mt);
              p > exp_mod_normal_large_p) {
            if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
              lp += log(-expm1(log_b) + g * (e_mt - e_mu));
            } else if (const auto lo = -t;
                       lo >= exp_mod_normal_erfcx_tail_start) {
              const auto hi = -u;
              const auto w0 = inv(lo * lo);
              const auto w1 = inv(hi * hi);
              T_partials a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
              T_partials pw0 = 1, pw1 = 1, sum_k = 0;
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
              lp += log(
                  -expm1(log_b)
                  + g
                        * (e_mt
                           * (d
                              * (1
                                 - lo
                                       * (SQRT_PI * (-w0 * ((hi + lo) * w1))
                                          * (dd_a * q0 - a0 * dd_q) / (q1 * q0)
                                          / (1
                                             + SQRT_PI * w0 * pq.first
                                                   / pq.second)))
                              / hi)));
            } else {
              const auto mid = -u - 0.5 * d;
              T_partials e_km1 = erfcx(mid);
              T_partials e_k = erfcx_derivative(mid, e_km1);
              std::array<T_partials, 4> odd;
              odd[0] = e_k;
              for (int k = 1; k < 7; ++k) {
                e_km1 = std::exchange(e_k, 2 * mid * e_k + 2 * k * e_km1);
                if (k % 2 == 0) {
                  odd[k / 2] = e_k;
                }
              }
              const auto h = 0.5 * d;
              const auto h2 = h * h;
              lp += log(
                  -expm1(log_b)
                  + g
                        * (d
                           * -(odd[0]
                               + h2
                                     * (odd[1] / 6
                                        + h2
                                              * (odd[2] / 120
                                                 + h2 * odd[3] / 5040)))));
            }
          } else {
            lp += log1p(-p);
          }
        }
      }
      return lp;
    }
  }
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
