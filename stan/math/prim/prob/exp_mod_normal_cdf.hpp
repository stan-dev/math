#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/exp_mod_normal_check_args.hpp>
#include <stan/math/prim/prob/exp_mod_normal_lcdf.hpp>
#include <array>
#include <cstddef>
#include <utility>

namespace stan {
namespace math {
namespace internal {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf_impl(
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
    return 1.0;
  } else if (unlikely(y_inf)) {
    if (any(y_val == NEGATIVE_INFTY)) {
      return 0.0;
    } else {
      return 1.0;
    }
  } else {
    scalar_seq_view<std::decay_t<decltype(y_val)>> y_vec(y_val);
    scalar_seq_view<std::decay_t<decltype(mu_val)>> mu_vec(mu_val);
    scalar_seq_view<std::decay_t<decltype(sigma_val)>> sigma_vec(sigma_val);
    scalar_seq_view<std::decay_t<decltype(lambda_val)>> lambda_vec(lambda_val);
    const std::size_t N = max_size(y_ref, mu_ref, sigma_ref, lambda_ref);
    T_partials cdf = 1.0;
    if constexpr (!any_autodiff) {
      for (std::size_t n = 0; n < N; ++n) {
        const auto sigma_n = sigma_vec[n];
        const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
        const auto v = lambda_vec[n] * sigma_n;
        const auto u = -z * INV_SQRT_TWO;
        const auto d = v * INV_SQRT_TWO;
        const auto t = u + d;
        const auto log_g = LOG_HALF - (0.5 * z) * z;
        if (z <= 0) {
          const auto e_t = erfcx(t);
          const auto e_u = erfcx(u);
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            cdf *= exp(log_g) * (e_u - e_t);
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            cdf *= exp(log_g)
                   * (e_u
                      * (d * (1 - u * exp_mod_normal_erfcx_tail_h_per_d(u, t))
                         / t));
          } else {
            const auto m = u + 0.5 * d;
            std::array<T_partials, 8> e;
            e[0] = erfcx(m);
            e[1] = erfcx_derivative(m, e[0]);
            for (std::size_t k = 1; k < 7; ++k) {
              e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                e[1] + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
            cdf *= exp(log_g) * (d * neg_s);
          }
        } else {
          const auto log_b = d * (2 * u + d);
          if (t < -6.5 && log_b < -0.2877) {
            cdf *= 1 - exp(log_b);
          } else {
            const auto g = exp(log_g);
            const auto e_mu = erfcx(-u);
            const auto q = g * e_mu;
            T_partials e_t = 0;
            T_partials e_mt = 0;
            T_partials b;
            if (t >= 0) {
              e_t = erfcx(t);
              b = g * e_t;
            } else {
              e_mt = erfcx(-t);
              b = exp(log_b) - g * e_mt;
            }
            const auto p = q + b;
            const bool small_f = p > 0.75;
            if (small_f && t >= 0) {
              const auto e_u = erfcx(u);
              if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
                cdf *= exp(log_g) * (e_u - e_t);
              } else if (u >= exp_mod_normal_erfcx_tail_start) {
                cdf *= exp(log_g)
                       * (e_u
                          * (d
                             * (1 - u * exp_mod_normal_erfcx_tail_h_per_d(u, t))
                             / t));
              } else {
                const auto m = u + 0.5 * d;
                std::array<T_partials, 8> e;
                e[0] = erfcx(m);
                e[1] = erfcx_derivative(m, e[0]);
                for (std::size_t k = 1; k < 7; ++k) {
                  e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
                }
                const auto h = 0.5 * d;
                const auto h2 = h * h;
                const auto neg_s = -(
                    e[1]
                    + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
                cdf *= exp(log_g) * (d * neg_s);
              }
            } else if (small_f) {
              const auto lo = -t;
              const auto hi = -u;
              if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
                cdf *= -expm1(log_b) + g * (e_mt - e_mu);
              } else if (lo >= exp_mod_normal_erfcx_tail_start) {
                cdf *= -expm1(log_b)
                       + g
                             * (e_mt
                                * (d
                                   * (1
                                      - lo
                                            * exp_mod_normal_erfcx_tail_h_per_d(
                                                lo, hi))
                                   / hi));
              } else {
                const auto m = hi - 0.5 * d;
                std::array<T_partials, 8> e;
                e[0] = erfcx(m);
                e[1] = erfcx_derivative(m, e[0]);
                for (std::size_t k = 1; k < 7; ++k) {
                  e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
                }
                const auto h = 0.5 * d;
                const auto h2 = h * h;
                const auto neg_s = -(
                    e[1]
                    + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
                cdf *= -expm1(log_b) + g * (d * neg_s);
              }
            } else {
              cdf *= 1 - p;
            }
          }
        }
      }
      return cdf;
    } else {
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
          const auto e_t = erfcx(t);
          const auto e_u = erfcx(u);
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            const auto delta = e_u - e_t;
            cdf *= exp(log_g) * delta;
            const auto scale = inv(delta);
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
            const auto delta
                = e_u
                  * (d * (1 - u * exp_mod_normal_erfcx_tail_h_per_d(u, t)) / t);
            cdf *= exp(log_g) * delta;
            const auto scale = inv(delta);
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
            const auto m = u + 0.5 * d;
            std::array<T_partials, 8> e;
            e[0] = erfcx(m);
            e[1] = erfcx_derivative(m, e[0]);
            for (std::size_t k = 1; k < 7; ++k) {
              e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                e[1] + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
            const auto delta = d * neg_s;
            cdf *= exp(log_g) * delta;
            const auto scale = inv(delta);
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
          const auto log_b = d * (2 * u + d);
          if (t < -6.5 && log_b < -0.2877) {
            const auto b = exp(log_b);
            const auto r = b / (1 - b);
            cdf *= 1 - b;
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
          } else {
            const auto g = exp(log_g);
            const auto e_mu = erfcx(-u);
            const auto q = g * e_mu;
            T_partials e_t = 0;
            T_partials e_mt = 0;
            T_partials b;
            if (t >= 0) {
              e_t = erfcx(t);
              b = g * e_t;
            } else {
              e_mt = erfcx(-t);
              b = exp(log_b) - g * e_mt;
            }
            const auto p = q + b;
            const bool small_f = p > 0.75;
            if (small_f && t >= 0) {
              const auto e_u = erfcx(u);
              if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
                const auto delta = e_u - e_t;
                cdf *= exp(log_g) * delta;
                const auto scale = inv(delta);
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
                const auto delta
                    = e_u
                      * (d * (1 - u * exp_mod_normal_erfcx_tail_h_per_d(u, t))
                         / t);
                cdf *= exp(log_g) * delta;
                const auto scale = inv(delta);
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
                const auto m = u + 0.5 * d;
                std::array<T_partials, 8> e;
                e[0] = erfcx(m);
                e[1] = erfcx_derivative(m, e[0]);
                for (std::size_t k = 1; k < 7; ++k) {
                  e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
                }
                const auto h = 0.5 * d;
                const auto h2 = h * h;
                const auto neg_s = -(
                    e[1]
                    + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
                const auto delta = d * neg_s;
                cdf *= exp(log_g) * delta;
                const auto scale = inv(delta);
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
            } else if (small_f) {
              const auto lo = -t;
              const auto hi = -u;
              if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
                const auto f = -expm1(log_b) + g * (e_mt - e_mu);
                cdf *= f;
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
              } else if (lo >= exp_mod_normal_erfcx_tail_start) {
                const auto f
                    = -expm1(log_b)
                      + g
                            * (e_mt
                               * (d
                                  * (1
                                     - lo
                                           * exp_mod_normal_erfcx_tail_h_per_d(
                                               lo, hi))
                                  / hi));
                cdf *= f;
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
                const auto m = hi - 0.5 * d;
                std::array<T_partials, 8> e;
                e[0] = erfcx(m);
                e[1] = erfcx_derivative(m, e[0]);
                for (std::size_t k = 1; k < 7; ++k) {
                  e[k + 1] = 2 * m * e[k] + 2 * k * e[k - 1];
                }
                const auto h = 0.5 * d;
                const auto h2 = h * h;
                const auto neg_s = -(
                    e[1]
                    + h2 * (e[3] / 6 + h2 * (e[5] / 120 + h2 * e[7] / 5040)));
                const auto f = -expm1(log_b) + g * (d * neg_s);
                cdf *= f;
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
              cdf *= f;
              if (t >= 0) {
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
              } else {
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
        }
      }
      if (unlikely(cdf == 0)) {
        return make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref)
            .build(0.0);
      } else {
        if constexpr (is_autodiff_v<T_y>) {
          for (std::size_t n = 0; n < stan::math::size(y_ref); ++n) {
            partials<0>(ops_partials)[n] *= cdf;
          }
        }
        if constexpr (is_autodiff_v<T_loc>) {
          for (std::size_t n = 0; n < stan::math::size(mu_ref); ++n) {
            partials<1>(ops_partials)[n] *= cdf;
          }
        }
        if constexpr (is_autodiff_v<T_scale>) {
          for (std::size_t n = 0; n < stan::math::size(sigma_ref); ++n) {
            partials<2>(ops_partials)[n] *= cdf;
          }
        }
        if constexpr (is_autodiff_v<T_inv_scale>) {
          for (std::size_t n = 0; n < stan::math::size(lambda_ref); ++n) {
            partials<3>(ops_partials)[n] *= cdf;
          }
        }
        return ops_partials.build(cdf);
      }
    }
  }
}

}  // namespace internal

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  return internal::exp_mod_normal_cdf_impl(
      "exp_mod_normal_cdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma), std::forward<T_inv_scale>(lambda));
}

}  // namespace math
}  // namespace stan
#endif
