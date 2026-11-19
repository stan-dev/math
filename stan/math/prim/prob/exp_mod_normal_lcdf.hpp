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

template <typename T>
inline auto exp_mod_normal_erfcx_h(T&& x) {
  const auto w = inv(x * x);
  const auto pq = erfcx_tail_polynomials<true>(w);
  return SQRT_PI * w * pq.first / pq.second;
}

template <typename T1, typename T2, typename T3>
inline auto exp_mod_normal_erfcx_h_diff(const T1& x, const T2& y,
                                        const T3& delta) {
  using T_return = return_type_t<T1, T2>;
  const auto w0 = inv(x * x);
  const auto w1 = inv(y * y);
  const auto dw = -(delta * w0) * ((y + x) * w1);
  T_return a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
  T_return pw0 = 1, pw1 = 1, sum_k = 0;
  for (std::size_t k = 0; k < erfcx_tail_p.size(); ++k) {
    q0 += erfcx_tail_q[k] * pw0;
    q1 += erfcx_tail_q[k] * pw1;
    if (k > 0) {
      dd_q += erfcx_tail_q[k] * sum_k;
    }
    sum_k = sum_k * w0 + pw1;
    a0 += erfcx_tail_p[k] * pw0 * w0;
    dd_a += erfcx_tail_p[k] * sum_k;
    pw0 *= w0;
    pw1 *= w1;
  }
  return SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
}

template <typename T>
struct exp_mod_normal_erfcx_diff_result {
  T delta;
  T log_delta;
  T d_e_hi_over_delta;
};

template <typename T>
inline T exp_mod_normal_log1p_over_x(const T& x) {
  return x == 0 ? T(1) : T(log1p(x) / x);
}

template <bool upper_exact, bool need_log, bool need_ratio, typename T1,
          typename T2, typename T3, typename T4, typename T5>
inline auto exp_mod_normal_erfcx_diff_close(const T1& lo, const T2& hi,
                                            const T3& d, const T4& e_lo,
                                            const T5& e_hi) {
  using T_return = return_type_t<T1, T2, T3, T4, T5>;
  if (lo >= exp_mod_normal_erfcx_tail_start) {
    const auto x_end = upper_exact ? T_return(-d / hi) : T_return(d / lo);
    const auto l_ratio_x
        = upper_exact ? T_return(-exp_mod_normal_log1p_over_x(x_end) / hi)
                      : T_return(-exp_mod_normal_log1p_over_x(x_end) / lo);
    const auto h_per_d = exp_mod_normal_erfcx_h_diff(lo, hi, 1.0)
                         / (1 + exp_mod_normal_erfcx_h(lo));
    const auto h_term = d * h_per_d;
    const auto l_total
        = l_ratio_x + h_per_d * exp_mod_normal_log1p_over_x(h_term);
    const auto log_r = d * l_total;
    const auto em = log_r == 0 ? T_return(1) : T_return(expm1(log_r) / log_r);
    const auto one_minus_r_over_d = -l_total * em;
    const auto one_minus_r = d * one_minus_r_over_d;
    return exp_mod_normal_erfcx_diff_result<T_return>{
        e_lo * one_minus_r,
        need_log ? T_return(log(e_lo) + log(d) + log(one_minus_r_over_d))
                 : T_return(0),
        need_ratio ? T_return(exp(log_r) / one_minus_r_over_d) : T_return(0)};
  }
  const auto m = upper_exact ? T_return(hi - 0.5 * d) : T_return(lo + 0.5 * d);
  T_return e_km1 = erfcx(m);
  T_return e_k = erfcx_derivative(m, e_km1);
  std::array<T_return, 4> odd;
  odd[0] = e_k;
  for (int k = 1; k < 7; ++k) {
    const auto e_kp1 = 2 * m * e_k + 2 * k * e_km1;
    e_km1 = e_k;
    e_k = e_kp1;
    if (k % 2 == 0) {
      odd[k / 2] = e_k;
    }
  }
  const auto h = 0.5 * d;
  const auto h2 = h * h;
  const auto neg_s = -(
      odd[0] + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
  return exp_mod_normal_erfcx_diff_result<T_return>{
      d * neg_s, need_log ? T_return(log(d) + log(neg_s)) : T_return(0),
      need_ratio ? T_return(e_hi / neg_s) : T_return(0)};
}

template <typename T1, typename T2>
inline bool exp_mod_normal_use_ratio(const T1& delta, const T2& hi) {
  return !(delta > 1e-300 && hi < 1e150);
}

template <typename T1, typename T2>
inline auto exp_mod_normal_erfcx_log_slope(const T1& t, const T2& e_t) {
  using T_return = return_type_t<T1, T2>;
  return t > 1e150 ? T_return(-inv(t))
                   : T_return(erfcx_derivative(t, e_t) / e_t);
}

template <bool upper_exact, bool need_log, bool need_ratio = false, typename T1,
          typename T2, typename T3, typename T4, typename T5>
inline auto exp_mod_normal_erfcx_diff(const T1& lo, const T2& hi, const T3& d,
                                      const T4& e_lo, const T5& e_hi) {
  using T_return = return_type_t<T1, T2, T3, T4, T5>;
  if (e_hi < exp_mod_normal_erfcx_close_ratio * e_lo) {
    const auto delta = e_lo - e_hi;
    return exp_mod_normal_erfcx_diff_result<T_return>{
        delta, need_log ? T_return(log(delta)) : T_return(0), 0};
  }
  return exp_mod_normal_erfcx_diff_close<upper_exact, need_log, need_ratio>(
      lo, hi, d, e_lo, e_hi);
}

template <typename T>
inline T exp_mod_normal_log1m(const T& x) {
  return x < 1e-17 ? T(-x) : T(log1p(-x));
}

template <typename T1, typename T2, typename T3, typename T4>
struct exp_mod_normal_lcdf_terms {
  T1 value;
  T2 d_z;
  T3 d_v;
  T4 d_s;
};

template <typename T1, typename T2, typename T3, typename T4>
exp_mod_normal_lcdf_terms(T1, T2, T3, T4)
    -> exp_mod_normal_lcdf_terms<T1, T2, T3, T4>;

template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
inline auto exp_mod_normal_lcdf_terms_from_ratio(const T1& value, const T2& u,
                                                 const T3& t, const T4& v,
                                                 const T5& e_t, const T6& dl) {
  using T_return = return_type_t<T2, T3, T4, T5>;
  const auto p = dl.d_e_hi_over_delta > 0
                     ? T_return(dl.d_e_hi_over_delta)
                     : T_return(v * INV_SQRT_TWO * (e_t / dl.delta));
  const auto s = exp_mod_normal_erfcx_log_slope(t, e_t);
  return exp_mod_normal_lcdf_terms{value, SQRT_TWO * p, -(s * p) / v,
                                   p * (2 * u - s)};
}

template <bool calc_grad, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6, typename T7, typename T8, typename T9>
inline auto exp_mod_normal_lccdf_tiny(const T1& z, const T2& v, const T3& u,
                                      const T4& t, const T5& log_g,
                                      const T6& log_b, const T7& e_t,
                                      const T8& e_mu, const T9& e_mt) {
  using T_return = return_type_t<T1, T2, T3, T4, T5, T6, T7, T8, T9>;
  const auto log_b_t = t >= 0
                           ? T_return(log_g + log(e_t))
                           : T_return(log_b + log1p(-0.5 * exp(-t * t) * e_mt));
  const auto log_p = log_sum_exp(T_return(log_g + log(e_mu)), log_b_t);
  if constexpr (!calc_grad) {
    return log_p;
  } else {
    const auto b_p = exp(log_b_t - log_p);
    const auto g_p = exp(log_g - log_p);
    if (t >= 0) {
      const auto de_t = erfcx_derivative(t, e_t);
      return exp_mod_normal_lcdf_terms{
          log_p, -v * b_p, g_p * de_t * INV_SQRT_TWO,
          -v * g_p * SQRT_TWO * (u * e_t - 0.5 * de_t)};
    }
    const auto phi_p = g_p * SQRT_TWO_OVER_SQRT_PI;
    return exp_mod_normal_lcdf_terms{log_p, -v * b_p, -(phi_p - (v - z) * b_p),
                                     -v * (phi_p - v * b_p)};
  }
}

template <typename T>
struct exp_mod_normal_middle_terms {
  T g;
  T e_mu;
  T e_t;
  T e_mt;
  T b;
  T p;
};

template <typename T1, typename T2, typename T3, typename T4>
inline auto exp_mod_normal_middle(const T1& u, const T2& t, const T3& log_g,
                                  const T4& log_b) {
  using T_return = return_type_t<T1, T2, T3, T4>;
  const auto g = exp(log_g);
  const auto e_mu = erfcx(-u);
  const auto q = g * e_mu;
  T_return e_t = 0;
  T_return e_mt = 0;
  T_return b;
  if (t >= 0) {
    e_t = erfcx(t);
    b = g * e_t;
  } else {
    e_mt = erfcx(-t);
    b = exp(log_b) - g * e_mt;
  }
  return exp_mod_normal_middle_terms<T_return>{g, e_mu, e_t, e_mt, b, q + b};
}

template <typename T>
struct exp_mod_normal_f_partials {
  T f_z;
  T f_v;
  T f_s;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline auto exp_mod_normal_partials_of_f(
    const T1& z, const T2& v, const T3& u, const T4& t,
    const exp_mod_normal_middle_terms<T5>& m) {
  using T_return = return_type_t<T1, T2, T3, T4, T5>;
  if (t >= 0) {
    const auto de_t = erfcx_derivative(t, m.e_t);
    return exp_mod_normal_f_partials<T_return>{
        v * m.b, -m.g * de_t * INV_SQRT_TWO,
        v * m.g * SQRT_TWO * (u * m.e_t - 0.5 * de_t)};
  }
  const auto phi = m.g * SQRT_TWO_OVER_SQRT_PI;
  return exp_mod_normal_f_partials<T_return>{v * m.b, phi - (v - z) * m.b,
                                             v * (phi - v * m.b)};
}

template <bool upper, bool calc_grad, typename T1, typename T2>
inline auto exp_mod_normal_lcdf_kernel(T1&& z, T2&& v) {
  if constexpr (upper) {
    if constexpr (calc_grad) {
      using T_return = return_type_t<T1, T2>;
      const auto u = -z * INV_SQRT_TWO;
      const auto d = v * INV_SQRT_TWO;
      const auto t = u + d;
      const auto log_g = LOG_HALF - (0.5 * z) * z;
      if (z <= 0) {
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        T_return delta;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          delta = e_u - e_t;
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = d / u;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / u;
          const auto w0 = inv(u * u);
          const auto w1 = inv(t * t);
          const auto dw = -w0 * ((t + u) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_u * one_minus_r;
        } else {
          const auto mid = u + 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        const auto g = exp(log_g);
        const auto f = g * delta;
        const auto value = log1p(-f);
        const auto scale = -g / (1 - f);
        const auto de_t = erfcx_derivative(t, e_t);
        return exp_mod_normal_lcdf_terms{
            value, scale * v * e_t, -scale * de_t * INV_SQRT_TWO,
            scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t)};
      }
      const auto log_b = d * (2 * u + d);
      if (t < exp_mod_normal_b_only_t && log_b < exp_mod_normal_log_large_p) {
        return exp_mod_normal_lcdf_terms{log_b, T_return(-v), v - z,
                                         T_return(v * v)};
      }
      const auto e_mu = erfcx(-u);
      if (t >= 0) {
        const auto g = exp(log_g);
        const auto e_t = erfcx(t);
        const auto e_sum = e_mu + e_t;
        const auto p = g * e_sum;
        if (p > exp_mod_normal_large_p) {
          const auto e_u = erfcx(u);
          T_return delta;
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            delta = e_u - e_t;
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto x_end = d / u;
            T_return log1p_over_x_end = 1;
            if (x_end != 0) {
              log1p_over_x_end = log1p(x_end) / x_end;
            }
            const auto l_ratio_x = -log1p_over_x_end / u;
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            const auto dw = -w0 * ((t + u) * w1);
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
            const auto h_diff_per_d
                = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
            const auto pq = erfcx_tail_polynomials<true>(w0);
            const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
            const auto h_per_d = h_diff_per_d / (1 + h_lo);
            const auto h_term = d * h_per_d;
            T_return log1p_over_h_term = 1;
            if (h_term != 0) {
              log1p_over_h_term = log1p(h_term) / h_term;
            }
            const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
            const auto log_r = d * l_total;
            T_return em = 1;
            if (log_r != 0) {
              em = expm1(log_r) / log_r;
            }
            const auto one_minus_r_over_d = -l_total * em;
            const auto one_minus_r = d * one_minus_r_over_d;
            delta = e_u * one_minus_r;
          } else {
            const auto mid = u + 0.5 * d;
            T_return e_km1 = erfcx(mid);
            T_return e_k = erfcx_derivative(mid, e_km1);
            std::array<T_return, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
              e_km1 = e_k;
              e_k = e_kp1;
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                odd[0]
                + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
            delta = d * neg_s;
          }
          const auto f = g * delta;
          const auto value = log1p(-f);
          const auto scale = -g / (1 - f);
          const auto de_t = erfcx_derivative(t, e_t);
          return exp_mod_normal_lcdf_terms{
              value, scale * v * e_t, -scale * de_t * INV_SQRT_TWO,
              scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t)};
        }
        const auto de_t = erfcx_derivative(t, e_t);
        return exp_mod_normal_lcdf_terms{
            log_g + log(e_sum), -v * e_t / e_sum, de_t * INV_SQRT_TWO / e_sum,
            -v * SQRT_TWO * (u * e_t - 0.5 * de_t) / e_sum};
      }
      const auto e_mt = erfcx(-t);
      const auto half_exp_mt2 = 0.5 * exp(-t * t);
      const auto rel_tail = half_exp_mt2 * (e_mt - e_mu);
      const auto p_over_b = 1 - rel_tail;
      const auto exp_log_b = exp(log_b);
      const auto p = exp_log_b * p_over_b;
      if (p > exp_mod_normal_large_p) {
        const auto lo = -t;
        const auto hi = -u;
        T_return delta;
        if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
          delta = e_mt - e_mu;
        } else if (lo >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = -d / hi;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / hi;
          const auto w0 = inv(lo * lo);
          const auto w1 = inv(hi * hi);
          const auto dw = -w0 * ((hi + lo) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_mt * one_minus_r;
        } else {
          const auto mid = hi - 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        const auto g = exp_log_b * half_exp_mt2;
        const auto b = exp_log_b * (1 - half_exp_mt2 * e_mt);
        const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
        const auto f = -expm1(log_b) + g * delta;
        const auto divisor = -(1 - f);
        return exp_mod_normal_lcdf_terms{log1p(-f), v * b / divisor,
                                         (phi - (v - z) * b) / divisor,
                                         v * (phi - v * b) / divisor};
      }
      const auto b_p = (1 - half_exp_mt2 * e_mt) / p_over_b;
      const auto phi_p = half_exp_mt2 * SQRT_TWO_OVER_SQRT_PI / p_over_b;
      return exp_mod_normal_lcdf_terms{log_b + log1p(-rel_tail), -v * b_p,
                                       -(phi_p - (v - z) * b_p),
                                       -v * (phi_p - v * b_p)};
    } else {
      using T_return = return_type_t<T1, T2>;
      const auto u = -z * INV_SQRT_TWO;
      const auto d = v * INV_SQRT_TWO;
      const auto t = u + d;
      const auto log_g = LOG_HALF - (0.5 * z) * z;
      if (z <= 0) {
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        T_return delta;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          delta = e_u - e_t;
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = d / u;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / u;
          const auto w0 = inv(u * u);
          const auto w1 = inv(t * t);
          const auto dw = -w0 * ((t + u) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_u * one_minus_r;
        } else {
          const auto mid = u + 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        const auto f = exp(log_g) * delta;
        return T_return(log1p(-f));
      }
      const auto log_b = d * (2 * u + d);
      if (t < exp_mod_normal_b_only_t && log_b < exp_mod_normal_log_large_p) {
        return log_b;
      }
      const auto e_mu = erfcx(-u);
      if (t >= 0) {
        const auto g = exp(log_g);
        const auto e_t = erfcx(t);
        const auto e_sum = e_mu + e_t;
        const auto p = g * e_sum;
        if (p > exp_mod_normal_large_p) {
          const auto e_u = erfcx(u);
          T_return delta;
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            delta = e_u - e_t;
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto x_end = d / u;
            T_return log1p_over_x_end = 1;
            if (x_end != 0) {
              log1p_over_x_end = log1p(x_end) / x_end;
            }
            const auto l_ratio_x = -log1p_over_x_end / u;
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            const auto dw = -w0 * ((t + u) * w1);
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
            const auto h_diff_per_d
                = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
            const auto pq = erfcx_tail_polynomials<true>(w0);
            const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
            const auto h_per_d = h_diff_per_d / (1 + h_lo);
            const auto h_term = d * h_per_d;
            T_return log1p_over_h_term = 1;
            if (h_term != 0) {
              log1p_over_h_term = log1p(h_term) / h_term;
            }
            const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
            const auto log_r = d * l_total;
            T_return em = 1;
            if (log_r != 0) {
              em = expm1(log_r) / log_r;
            }
            const auto one_minus_r_over_d = -l_total * em;
            const auto one_minus_r = d * one_minus_r_over_d;
            delta = e_u * one_minus_r;
          } else {
            const auto mid = u + 0.5 * d;
            T_return e_km1 = erfcx(mid);
            T_return e_k = erfcx_derivative(mid, e_km1);
            std::array<T_return, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
              e_km1 = e_k;
              e_k = e_kp1;
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                odd[0]
                + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
            delta = d * neg_s;
          }
          const auto f = g * delta;
          return T_return(log1p(-f));
        }
        return log_g + log(e_sum);
      }
      const auto e_mt = erfcx(-t);
      const auto half_exp_mt2 = 0.5 * exp(-t * t);
      const auto rel_tail = half_exp_mt2 * (e_mt - e_mu);
      const auto exp_log_b = exp(log_b);
      const auto p = exp_log_b * (1 - rel_tail);
      if (p > exp_mod_normal_large_p) {
        const auto lo = -t;
        const auto hi = -u;
        T_return delta;
        if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
          delta = e_mt - e_mu;
        } else if (lo >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = -d / hi;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / hi;
          const auto w0 = inv(lo * lo);
          const auto w1 = inv(hi * hi);
          const auto dw = -w0 * ((hi + lo) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_mt * one_minus_r;
        } else {
          const auto mid = hi - 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        const auto g = exp_log_b * half_exp_mt2;
        const auto f = -expm1(log_b) + g * delta;
        return log1p(-f);
      }
      return log_b + log1p(-rel_tail);
    }
  } else {
    if constexpr (calc_grad) {
      using T_return = return_type_t<T1, T2>;
      const auto u = -z * INV_SQRT_TWO;
      const auto d = v * INV_SQRT_TWO;
      const auto t = u + d;
      const auto log_g = LOG_HALF - (0.5 * z) * z;
      if (z <= 0) {
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        T_return log_delta;
        T_return d_e_t_over_delta;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          const auto delta = e_u - e_t;
          log_delta = log(delta);
          d_e_t_over_delta = d * (e_t / delta);
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = d / u;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / u;
          const auto w0 = inv(u * u);
          const auto w1 = inv(t * t);
          const auto dw = -w0 * ((t + u) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          log_delta = log(e_u) + log(d) + log(one_minus_r_over_d);
          d_e_t_over_delta = exp(log_r) / one_minus_r_over_d;
        } else {
          const auto mid = u + 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          log_delta = log(d) + log(neg_s);
          d_e_t_over_delta = e_t / neg_s;
        }
        T_return s;
        if (t > 1e150) {
          s = -inv(t);
        } else {
          s = erfcx_derivative(t, e_t) / e_t;
        }
        return exp_mod_normal_lcdf_terms{
            log_g + log_delta, SQRT_TWO * d_e_t_over_delta,
            -(s * d_e_t_over_delta) / v, d_e_t_over_delta * (2 * u - s)};
      }
      const auto log_b = d * (2 * u + d);
      if (t < exp_mod_normal_b_only_t && log_b < exp_mod_normal_log_large_p) {
        const auto b = exp(log_b);
        const auto r = b / (1 - b);
        return exp_mod_normal_lcdf_terms{log1p(-b), v * r, (z - v) * r,
                                         -v * v * r};
      }
      const auto g = exp(log_g);
      const auto e_mu = erfcx(-u);
      const auto q = g * e_mu;
      if (t >= 0) {
        const auto e_t = erfcx(t);
        const auto b = g * e_t;
        const auto p = q + b;
        if (p > exp_mod_normal_large_p) {
          const auto e_u = erfcx(u);
          T_return log_delta;
          T_return d_e_t_over_delta;
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            const auto delta = e_u - e_t;
            log_delta = log(delta);
            d_e_t_over_delta = d * (e_t / delta);
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto x_end = d / u;
            T_return log1p_over_x_end = 1;
            if (x_end != 0) {
              log1p_over_x_end = log1p(x_end) / x_end;
            }
            const auto l_ratio_x = -log1p_over_x_end / u;
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            const auto dw = -w0 * ((t + u) * w1);
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
            const auto h_diff_per_d
                = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
            const auto pq = erfcx_tail_polynomials<true>(w0);
            const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
            const auto h_per_d = h_diff_per_d / (1 + h_lo);
            const auto h_term = d * h_per_d;
            T_return log1p_over_h_term = 1;
            if (h_term != 0) {
              log1p_over_h_term = log1p(h_term) / h_term;
            }
            const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
            const auto log_r = d * l_total;
            T_return em = 1;
            if (log_r != 0) {
              em = expm1(log_r) / log_r;
            }
            const auto one_minus_r_over_d = -l_total * em;
            log_delta = log(e_u) + log(d) + log(one_minus_r_over_d);
            d_e_t_over_delta = exp(log_r) / one_minus_r_over_d;
          } else {
            const auto mid = u + 0.5 * d;
            T_return e_km1 = erfcx(mid);
            T_return e_k = erfcx_derivative(mid, e_km1);
            std::array<T_return, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
              e_km1 = e_k;
              e_k = e_kp1;
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                odd[0]
                + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
            log_delta = log(d) + log(neg_s);
            d_e_t_over_delta = e_t / neg_s;
          }
          T_return s;
          if (t > 1e150) {
            s = -inv(t);
          } else {
            s = erfcx_derivative(t, e_t) / e_t;
          }
          return exp_mod_normal_lcdf_terms{
              log_g + log_delta, SQRT_TWO * d_e_t_over_delta,
              -(s * d_e_t_over_delta) / v, d_e_t_over_delta * (2 * u - s)};
        }
        const auto de_t = erfcx_derivative(t, e_t);
        const auto f_z = v * b;
        const auto f_v = -g * de_t * INV_SQRT_TWO;
        const auto f_s = v * g * SQRT_TWO * (u * e_t - 0.5 * de_t);
        const auto f = 1 - p;
        return exp_mod_normal_lcdf_terms{log1p(-p), f_z / f, f_v / f, f_s / f};
      }
      const auto e_mt = erfcx(-t);
      const auto b = exp(log_b) - g * e_mt;
      const auto p = q + b;
      const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
      const auto f_z = v * b;
      const auto f_v = phi - (v - z) * b;
      const auto f_s = v * (phi - v * b);
      if (p > exp_mod_normal_large_p) {
        const auto lo = -t;
        const auto hi = -u;
        T_return delta;
        if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
          delta = e_mt - e_mu;
        } else if (lo >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = -d / hi;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / hi;
          const auto w0 = inv(lo * lo);
          const auto w1 = inv(hi * hi);
          const auto dw = -w0 * ((hi + lo) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_mt * one_minus_r;
        } else {
          const auto mid = hi - 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        const auto f = -expm1(log_b) + g * delta;
        return exp_mod_normal_lcdf_terms{log(f), f_z / f, f_v / f, f_s / f};
      }
      const auto f = 1 - p;
      return exp_mod_normal_lcdf_terms{log1p(-p), f_z / f, f_v / f, f_s / f};
    } else {
      using T_return = return_type_t<T1, T2>;
      const auto u = -z * INV_SQRT_TWO;
      const auto d = v * INV_SQRT_TWO;
      const auto t = u + d;
      const auto log_g = LOG_HALF - (0.5 * z) * z;
      if (z <= 0) {
        const auto e_t = erfcx(t);
        const auto e_u = erfcx(u);
        T_return log_delta;
        if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
          log_delta = log(e_u - e_t);
        } else if (u >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = d / u;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / u;
          const auto w0 = inv(u * u);
          const auto w1 = inv(t * t);
          const auto dw = -w0 * ((t + u) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          log_delta = log(e_u) + log(d) + log(one_minus_r_over_d);
        } else {
          const auto mid = u + 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          log_delta = log(d) + log(neg_s);
        }
        return log_g + log_delta;
      }
      const auto log_b = d * (2 * u + d);
      if (t < exp_mod_normal_b_only_t && log_b < exp_mod_normal_log_large_p) {
        const auto b = exp(log_b);
        return T_return(log1p(-b));
      }
      const auto g = exp(log_g);
      const auto e_mu = erfcx(-u);
      const auto q = g * e_mu;
      if (t >= 0) {
        const auto e_t = erfcx(t);
        const auto b = g * e_t;
        const auto p = q + b;
        if (p > exp_mod_normal_large_p) {
          const auto e_u = erfcx(u);
          T_return log_delta;
          if (e_t < exp_mod_normal_erfcx_close_ratio * e_u) {
            log_delta = log(e_u - e_t);
          } else if (u >= exp_mod_normal_erfcx_tail_start) {
            const auto x_end = d / u;
            T_return log1p_over_x_end = 1;
            if (x_end != 0) {
              log1p_over_x_end = log1p(x_end) / x_end;
            }
            const auto l_ratio_x = -log1p_over_x_end / u;
            const auto w0 = inv(u * u);
            const auto w1 = inv(t * t);
            const auto dw = -w0 * ((t + u) * w1);
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
            const auto h_diff_per_d
                = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
            const auto pq = erfcx_tail_polynomials<true>(w0);
            const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
            const auto h_per_d = h_diff_per_d / (1 + h_lo);
            const auto h_term = d * h_per_d;
            T_return log1p_over_h_term = 1;
            if (h_term != 0) {
              log1p_over_h_term = log1p(h_term) / h_term;
            }
            const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
            const auto log_r = d * l_total;
            T_return em = 1;
            if (log_r != 0) {
              em = expm1(log_r) / log_r;
            }
            const auto one_minus_r_over_d = -l_total * em;
            log_delta = log(e_u) + log(d) + log(one_minus_r_over_d);
          } else {
            const auto mid = u + 0.5 * d;
            T_return e_km1 = erfcx(mid);
            T_return e_k = erfcx_derivative(mid, e_km1);
            std::array<T_return, 4> odd;
            odd[0] = e_k;
            for (int k = 1; k < 7; ++k) {
              const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
              e_km1 = e_k;
              e_k = e_kp1;
              if (k % 2 == 0) {
                odd[k / 2] = e_k;
              }
            }
            const auto h = 0.5 * d;
            const auto h2 = h * h;
            const auto neg_s = -(
                odd[0]
                + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
            log_delta = log(d) + log(neg_s);
          }
          return log_g + log_delta;
        }
        return T_return(log1p(-p));
      }
      const auto e_mt = erfcx(-t);
      const auto b = exp(log_b) - g * e_mt;
      const auto p = q + b;
      if (p > exp_mod_normal_large_p) {
        const auto lo = -t;
        const auto hi = -u;
        T_return delta;
        if (e_mu < exp_mod_normal_erfcx_close_ratio * e_mt) {
          delta = e_mt - e_mu;
        } else if (lo >= exp_mod_normal_erfcx_tail_start) {
          const auto x_end = -d / hi;
          T_return log1p_over_x_end = 1;
          if (x_end != 0) {
            log1p_over_x_end = log1p(x_end) / x_end;
          }
          const auto l_ratio_x = -log1p_over_x_end / hi;
          const auto w0 = inv(lo * lo);
          const auto w1 = inv(hi * hi);
          const auto dw = -w0 * ((hi + lo) * w1);
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
          const auto h_diff_per_d
              = SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
          const auto pq = erfcx_tail_polynomials<true>(w0);
          const auto h_lo = SQRT_PI * w0 * pq.first / pq.second;
          const auto h_per_d = h_diff_per_d / (1 + h_lo);
          const auto h_term = d * h_per_d;
          T_return log1p_over_h_term = 1;
          if (h_term != 0) {
            log1p_over_h_term = log1p(h_term) / h_term;
          }
          const auto l_total = l_ratio_x + h_per_d * log1p_over_h_term;
          const auto log_r = d * l_total;
          T_return em = 1;
          if (log_r != 0) {
            em = expm1(log_r) / log_r;
          }
          const auto one_minus_r_over_d = -l_total * em;
          const auto one_minus_r = d * one_minus_r_over_d;
          delta = e_mt * one_minus_r;
        } else {
          const auto mid = hi - 0.5 * d;
          T_return e_km1 = erfcx(mid);
          T_return e_k = erfcx_derivative(mid, e_km1);
          std::array<T_return, 4> odd;
          odd[0] = e_k;
          for (int k = 1; k < 7; ++k) {
            const auto e_kp1 = 2 * mid * e_k + 2 * k * e_km1;
            e_km1 = e_k;
            e_k = e_kp1;
            if (k % 2 == 0) {
              odd[k / 2] = e_k;
            }
          }
          const auto h = 0.5 * d;
          const auto h2 = h * h;
          const auto neg_s = -(
              odd[0]
              + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
          delta = d * neg_s;
        }
        return log(-expm1(log_b) + g * delta);
      }
      return T_return(log1p(-p));
    }
  }
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
  if constexpr (!any_autodiff) {
    for (std::size_t n = 0; n < N; ++n) {
      const auto sigma_n = sigma_vec[n];
      const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
      const auto v = lambda_vec[n] * sigma_n;
      lp += exp_mod_normal_lcdf_kernel<upper, false>(z, v);
    }
    return lp;
  } else {
    auto ops_partials
        = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);
    for (std::size_t n = 0; n < N; ++n) {
      const auto sigma_n = sigma_vec[n];
      const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
      const auto v = lambda_vec[n] * sigma_n;
      const auto r = exp_mod_normal_lcdf_kernel<upper, true>(z, v);
      lp += r.value;
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += r.d_z / sigma_n;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= r.d_z / sigma_n;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n] += r.d_s / sigma_n;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<3>(ops_partials)[n] += r.d_v * sigma_n;
      }
    }
    return ops_partials.build(lp);
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
