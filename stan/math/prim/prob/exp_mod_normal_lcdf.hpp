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
#include <cstddef>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * For x >= 4, erfcx(x) = (1 + h(x)) / (sqrt(pi) x) with the Cody tail
 * rational: h = sqrt(pi) w P(w) / Q(w), w = 1 / x^2.
 */
template <typename T>
inline T exp_mod_normal_erfcx_h(const T& x) {
  const T w = inv(x * x);
  const auto pq = erfcx_tail_polynomials<erfcx_tail_p.size()>(
      erfcx_tail_p.crbegin(), erfcx_tail_q.crbegin(), w);
  return SQRT_PI * w * pq.first / pq.second;
}

/**
 * h(y) - h(x) for x, y >= 4, with y - x = delta given exactly. The difference
 * of the rational A(w) / Q(w), A(w) = w P(w), is written as
 * (w1 - w0) [DD_A Q(w0) - A(w0) DD_Q] / (Q(w1) Q(w0)) with the divided
 * differences DD of the polynomials, so nothing cancels.
 */
template <typename T>
inline T exp_mod_normal_erfcx_h_diff(const T& x, const T& y, const T& delta) {
  const T w0 = inv(x * x);
  const T w1 = inv(y * y);
  // Grouped so that no factor overflows: delta * (y + x) overflows for
  // x above about 1e154, where w0 and w1 are 0.
  const T dw = -(delta * w0) * ((y + x) * w1);
  // Powers w1^j w0^(k-1-j) summed over j, for k = 1..6.
  T a0 = 0, q0 = 0, q1 = 0, dd_a = 0, dd_q = 0;
  T pw0 = 1, pw1 = 1, sum_k = 0;  // sum_k = sum_{j<k} w1^j w0^(k-1-j)
  for (std::size_t k = 0; k < erfcx_tail_p.size(); ++k) {
    // A has coefficient p[k] at power k + 1; Q has q[k] at power k.
    q0 += erfcx_tail_q[k] * pw0;
    q1 += erfcx_tail_q[k] * pw1;
    if (k > 0) {
      dd_q += erfcx_tail_q[k] * sum_k;
    }
    sum_k = sum_k * w0 + pw1;  // now the sum for power k + 1
    a0 += erfcx_tail_p[k] * pw0 * w0;
    dd_a += erfcx_tail_p[k] * sum_k;
    pw0 *= w0;
    pw1 *= w1;
  }
  return SQRT_PI * dw * (dd_a * q0 - a0 * dd_q) / (q1 * q0);
}

/**
 * erfcx(lo) - erfcx(hi), and its log when need_log is set, for hi - lo = d > 0
 * with d given exactly; e_lo and e_hi are erfcx(lo) and erfcx(hi). The end
 * that is not computed from the other one is named by upper_exact; the
 * difference never uses a rounded hi - lo.
 *  - erfcx(hi) / erfcx(lo) < 0.99: the direct difference;
 *  - else for lo >= 4: the Cody tail form,
 *    erfcx(hi) / erfcx(lo) = (lo / hi) (1 + h(hi)) / (1 + h(lo));
 *  - otherwise: the midpoint series in the odd derivatives, from
 *    E^(k+1) = 2 m E^(k) + 2 k E^(k-1).
 */
template <typename T>
struct exp_mod_normal_erfcx_diff_result {
  T delta;
  T log_delta;
  // d erfcx(hi) / delta from the close cases when need_ratio is set, else 0:
  // finite where delta underflows and 1 / delta overflows.
  T d_e_hi_over_delta;
};

/** log1p(x) / x, 1 at x == 0. */
template <typename T>
inline T exp_mod_normal_log1p_over_x(const T& x) {
  return x == 0 ? T(1) : T(log1p(x) / x);
}

/** The cases with erfcx(hi) / erfcx(lo) >= 0.99; rare unless d is small. */
template <bool upper_exact, bool need_log, bool need_ratio, typename T>
STAN_COLD_PATH exp_mod_normal_erfcx_diff_result<T>
exp_mod_normal_erfcx_diff_close(const T& lo, const T& hi, const T& d,
                                const T& e_lo, const T& e_hi) {
  if (lo >= 4) {
    // log r = log(erfcx(hi) / erfcx(lo)) = d L, with L formed from
    // log1p(x) / x factors so that it stays finite when d / lo underflows.
    const T x_end = upper_exact ? T(-d / hi) : T(d / lo);
    const T l_ratio_x = upper_exact
                            ? T(-exp_mod_normal_log1p_over_x(x_end) / hi)
                            : T(-exp_mod_normal_log1p_over_x(x_end) / lo);
    const T h_per_d = exp_mod_normal_erfcx_h_diff(lo, hi, T(1.0))
                      / (1 + exp_mod_normal_erfcx_h(lo));
    // h_diff is linear in its delta argument.
    const T h_term = d * h_per_d;
    const T l_total = l_ratio_x + h_per_d * exp_mod_normal_log1p_over_x(h_term);
    const T log_r = d * l_total;
    // (1 - r) / d = -L expm1(log r) / log r
    const T em = log_r == 0 ? T(1) : T(expm1(log_r) / log_r);
    const T one_minus_r_over_d = -l_total * em;
    const T one_minus_r = d * one_minus_r_over_d;
    exp_mod_normal_erfcx_diff_result<T> res{e_lo * one_minus_r, 0, 0};
    if constexpr (need_log) {
      // Separate logs: the product can underflow for a very small d.
      res.log_delta = log(e_lo) + log(d) + log(one_minus_r_over_d);
    }
    if constexpr (need_ratio) {
      // d e_hi / (e_lo (1 - r)) = r / ((1 - r) / d)
      res.d_e_hi_over_delta = exp(log_r) / one_minus_r_over_d;
    }
    return res;
  }
  const T m = upper_exact ? T(hi - 0.5 * d) : T(lo + 0.5 * d);
  T e_km1 = erfcx(m);
  T e_k = erfcx_derivative(m, e_km1);
  // odd[j] = E^(2j+1)(m)
  T odd[4];
  odd[0] = e_k;
  for (int k = 1; k < 7; ++k) {
    const T e_kp1 = 2 * m * e_k + 2 * k * e_km1;
    e_km1 = e_k;
    e_k = e_kp1;
    if (k % 2 == 0) {
      odd[k / 2] = e_k;
    }
  }
  const T h = 0.5 * d;
  const T h2 = h * h;
  // delta = -d s with s > 0 for every d
  const T neg_s = -(
      odd[0] + h2 * (odd[1] / 6 + h2 * (odd[2] / 120 + h2 * odd[3] / 5040)));
  exp_mod_normal_erfcx_diff_result<T> res{d * neg_s, 0, 0};
  if constexpr (need_log) {
    res.log_delta = log(d) + log(neg_s);
  }
  if constexpr (need_ratio) {
    res.d_e_hi_over_delta = e_hi / neg_s;
  }
  return res;
}

/**
 * Whether the partials of log F need the ratio d erfcx(hi) / delta: 1 / delta
 * overflows below delta = 1e-300 (deep lower tail, or a tiny v), and
 * erfcx'(hi) underflows above hi = 1e150.
 */
template <typename T>
inline bool exp_mod_normal_use_ratio(const T& delta, const T& hi) {
  return !(delta > 1e-300 && hi < 1e150);
}

/**
 * erfcx'(t) / erfcx(t). Above t = 1e150, where erfcx' underflows, the ratio
 * is -1 / t to working precision.
 */
template <typename T>
inline T exp_mod_normal_erfcx_log_slope(const T& t, const T& e_t) {
  return t > 1e150 ? T(-inv(t)) : T(erfcx_derivative(t, e_t) / e_t);
}

template <bool upper_exact, bool need_log, bool need_ratio = false, typename T>
inline exp_mod_normal_erfcx_diff_result<T> exp_mod_normal_erfcx_diff(
    const T& lo, const T& hi, const T& d, const T& e_lo, const T& e_hi) {
  if (e_hi < 0.99 * e_lo) {
    const T delta = e_lo - e_hi;
    exp_mod_normal_erfcx_diff_result<T> res{delta, 0, 0};
    if constexpr (need_log) {
      res.log_delta = log(delta);
    }
    // d_e_hi_over_delta stays 0: here delta >= 0.01 erfcx(lo), and the
    // partials form the ratio from delta where they need it.
    return res;
  }
  return exp_mod_normal_erfcx_diff_close<upper_exact, need_log, need_ratio>(
      lo, hi, d, e_lo, e_hi);
}

/**
 * log(1 - x) for 0 <= x < 1. Below 1e-17 the result is -x, as the next term
 * x^2 / 2 is below the rounding error.
 */
template <typename T>
inline T exp_mod_normal_log1m(const T& x) {
  return x < 1e-17 ? T(-x) : T(log1p(-x));
}

/**
 * The log cdf (or the log ccdf when upper is set) of the exponentially
 * modified normal and its partials with respect to z, v and sigma d/dsigma at
 * fixed y, mu and lambda; z = (y - mu) / sigma, v = lambda sigma.
 */
template <typename T>
struct exp_mod_normal_lcdf_terms {
  T value;
  T d_z;
  T d_v;
  T d_s;
};

/**
 * The terms from F = g Delta, Delta = E(u) - E(t) with t >= 0, for F <= 1/2.
 * The partials of F are g times (v E(t), -E'(t) / sqrt 2,
 * v sqrt 2 (u E(t) - E'(t) / 2)). The log cdf divides them by F, the log ccdf
 * by -(1 - F). The log ccdf needs no log of Delta, and the cdf (log_value
 * unset) takes F = g Delta itself.
 * Where 1 / Delta overflows or E'(t) underflows (exp_mod_normal_use_ratio),
 * the partials of log F come from p = d E(t) / Delta and s = E'(t) / E(t)
 * with d = v / sqrt 2: (sqrt 2 p, -s p / v, p (2 u - s)).
 */
/**
 * The partials of log F (value field unset) from p = d E(t) / Delta and
 * s = E'(t) / E(t), for where exp_mod_normal_use_ratio holds.
 */
template <typename T>
STAN_COLD_PATH exp_mod_normal_lcdf_terms<T>
exp_mod_normal_lcdf_partials_from_ratio(
    const T& u, const T& t, const T& v, const T& e_t,
    const exp_mod_normal_erfcx_diff_result<T>& dl) {
  // The direct difference leaves the ratio 0; there delta is not small
  // against erfcx(u), so d e_t / delta has no overflow.
  const T p = dl.d_e_hi_over_delta > 0 ? dl.d_e_hi_over_delta
                                       : T(v * INV_SQRT_TWO * (e_t / dl.delta));
  const T s = exp_mod_normal_erfcx_log_slope(t, e_t);
  return {0, SQRT_TWO * p, -(s * p) / v, p * (2 * u - s)};
}

template <bool upper, bool calc_grad, bool log_value, typename T>
inline exp_mod_normal_lcdf_terms<T> exp_mod_normal_lcdf_from_delta(
    const T& u, const T& t, const T& v, const T& log_g, const T& e_t,
    const exp_mod_normal_erfcx_diff_result<T>& dl) {
  exp_mod_normal_lcdf_terms<T> r{0, 0, 0, 0};
  if constexpr (upper) {
    const T g = exp(log_g);
    const T f = g * dl.delta;
    r.value = exp_mod_normal_log1m(f);
    if constexpr (calc_grad) {
      const T scale = -g / (1 - f);
      const T de_t = erfcx_derivative(t, e_t);
      r.d_z = scale * v * e_t;
      r.d_v = -scale * de_t * INV_SQRT_TWO;
      r.d_s = scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t);
    }
    return r;
  }
  if constexpr (log_value) {
    r.value = log_g + dl.log_delta;
  } else {
    r.value = exp(log_g) * dl.delta;
  }
  if constexpr (calc_grad) {
    const T scale = inv(dl.delta);
    const T de_t = erfcx_derivative(t, e_t);
    r.d_z = scale * v * e_t;
    r.d_v = -scale * de_t * INV_SQRT_TWO;
    r.d_s = scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t);
    if (unlikely(exp_mod_normal_use_ratio(dl.delta, t))) {
      const auto q = exp_mod_normal_lcdf_partials_from_ratio(u, t, v, e_t, dl);
      r.d_z = q.d_z;
      r.d_v = q.d_v;
      r.d_s = q.d_s;
    }
  }
  return r;
}

/**
 * The log ccdf where P = Phi(-z) + b is near underflow (z near 38 and above),
 * from the logs of the terms. For t < 0, b = exp(d (2u + d)) (1 - Phi(v - z))
 * with Phi(v - z) < 1/2.
 */
template <bool calc_grad, typename T>
STAN_COLD_PATH exp_mod_normal_lcdf_terms<T> exp_mod_normal_lccdf_tiny(
    const T& z, const T& v, const T& u, const T& t, const T& log_g,
    const T& log_b, const T& e_t, const T& e_mu, const T& e_mt) {
  const T log_b_t = t >= 0 ? T(log_g + log(e_t))
                           : T(log_b + log1p(-0.5 * exp(-t * t) * e_mt));
  const T log_p = log_sum_exp(T(log_g + log(e_mu)), log_b_t);
  if constexpr (calc_grad) {
    const T b_p = exp(log_b_t - log_p);
    const T g_p = exp(log_g - log_p);
    if (t >= 0) {
      const T de_t = erfcx_derivative(t, e_t);
      return {log_p, -v * b_p, g_p * de_t * INV_SQRT_TWO,
              -v * g_p * SQRT_TWO * (u * e_t - 0.5 * de_t)};
    }
    const T phi_p = g_p * SQRT_TWO_OVER_SQRT_PI;
    return {log_p, -v * b_p, -(phi_p - (v - z) * b_p), -v * (phi_p - v * b_p)};
  } else {
    return {log_p, 0, 0, 0};
  }
}

/**
 * With u = -z / sqrt 2, d = v / sqrt 2, t = u + d, E = erfcx and
 * g = exp(-z^2 / 2) / 2:
 *   Phi(z) = g E(u),  b = exp(v^2 / 2 - v z) Phi(z - v) = g E(t),
 *   F = Phi(z) - b = g (E(u) - E(t)),  P = 1 - F = Phi(-z) + b,
 *   F_z = v b,  F_v = phi(z) - (v - z) b = -g E'(t) / sqrt 2,
 *   sigma F_sigma = -z F_z + v F_v = v (phi(z) - v b)
 *                 = v g sqrt 2 (u E(t) - E'(t) / 2).
 * The z^2 / 2 terms cancel analytically, so the tails need no difference of
 * large logs. For t < 0, b = exp(d (2u + d)) - g E(-t), as
 * E(t) = 2 exp(t^2) - E(-t). For z <= 0 (then F <= 1/2), and for z > 0 with
 * F < 1/4: F = g (E(u) - E(t)), or for t < 0,
 * F = -expm1(d (2u + d)) + g (E(-t) - E(-u)) with both terms positive.
 * Otherwise F = 1 - P with P = Phi(-z) + b.
 * With log_value unset (the cdf, lower tail only) the value is F itself and
 * the partials are still those of log F.
 */
template <bool upper, bool calc_grad, bool log_value, typename T>
inline exp_mod_normal_lcdf_terms<T> exp_mod_normal_lcdf_kernel(const T& z,
                                                               const T& v) {
  static_assert(log_value || !upper, "the ccdf has no plain value kernel");
  constexpr bool need_log = !upper && log_value;
  constexpr bool need_ratio = calc_grad && !upper;
  const T u = -z * INV_SQRT_TWO;
  const T d = v * INV_SQRT_TWO;
  const T t = u + d;
  const T log_g = LOG_HALF - (0.5 * z) * z;
  if (z <= 0) {
    // F <= Phi(z) <= 1/2.
    const T e_t = erfcx(t);
    const auto dl = exp_mod_normal_erfcx_diff<false, need_log, need_ratio>(
        u, t, d, T(erfcx(u)), e_t);
    return exp_mod_normal_lcdf_from_delta<upper, calc_grad, log_value>(
        u, t, v, log_g, e_t, dl);
  }
  const T log_b = d * (2 * u + d);
  if (t < -6.5 && log_b < -0.2877) {
    // Phi(-z), Phi(v - z) and phi(z) are below 3e-19 of b = exp(d (2u + d)),
    // and b < 3/4: P = b, and the partials of F are b (v, z - v, -v^2).
    if constexpr (upper) {
      return {log_b, -v, v - z, v * v};
    } else {
      const T b = exp(log_b);
      const T r = b / (1 - b);
      return {log_value ? exp_mod_normal_log1m(b) : T(1 - b), v * r,
              (z - v) * r, -v * v * r};
    }
  }
  const T g = exp(log_g);
  const T e_mu = erfcx(-u);
  const T q = g * e_mu;
  T e_t = 0;
  T e_mt = 0;
  T b;
  if (t >= 0) {
    e_t = erfcx(t);
    b = g * e_t;
  } else {
    e_mt = erfcx(-t);
    b = exp(log_b) - g * e_mt;
  }
  const T p = q + b;
  // F < 1/4: 1 - p would lose digits.
  const bool small_f = p > 0.75;
  if (small_f && t >= 0) {
    const auto dl = exp_mod_normal_erfcx_diff<false, need_log, need_ratio>(
        u, t, d, T(erfcx(u)), e_t);
    return exp_mod_normal_lcdf_from_delta<upper, calc_grad, log_value>(
        u, t, v, log_g, e_t, dl);
  }
  // The partials of F: f_z, f_v, f_s.
  T f_z = 0;
  T f_v = 0;
  T f_s = 0;
  if constexpr (calc_grad) {
    f_z = v * b;
    if (t >= 0) {
      const T de_t = erfcx_derivative(t, e_t);
      f_v = -g * de_t * INV_SQRT_TWO;
      f_s = v * g * SQRT_TWO * (u * e_t - 0.5 * de_t);
    } else {
      const T phi = g * SQRT_TWO_OVER_SQRT_PI;
      f_v = phi - (v - z) * b;
      f_s = v * (phi - v * b);
    }
  }
  if (small_f) {
    // t < 0: F = -expm1(d (2u + d)) + g (E(-t) - E(-u)).
    const auto dl
        = exp_mod_normal_erfcx_diff<true, false>(T(-t), T(-u), d, e_mt, e_mu);
    const T f = -expm1(log_b) + g * dl.delta;
    if constexpr (!upper) {
      return {log_value ? T(log(f)) : f, f_z / f, f_v / f, f_s / f};
    } else {
      const T one_m_f = 1 - f;
      return {log1p(-f), -f_z / one_m_f, -f_v / one_m_f, -f_s / one_m_f};
    }
  }
  // F >= 1/4: F = 1 - p loses at most 2 bits.
  if constexpr (!upper) {
    // log1p only where log F is small: for p >= 1/4, |log F| >= 0.29.
    const T f = 1 - p;
    if constexpr (!log_value) {
      return {f, f_z / f, f_v / f, f_s / f};
    }
    return {p < 0.25 ? exp_mod_normal_log1m(p) : T(log(f)), f_z / f, f_v / f,
            f_s / f};
  } else {
    if (p > 1e-300) {
      return {log(p), -f_z / p, -f_v / p, -f_s / p};
    }
    return exp_mod_normal_lccdf_tiny<calc_grad>(z, v, u, t, log_g, log_b, e_t,
                                                e_mu, e_mt);
  }
}

/**
 * log F, or log(1 - F) when upper is set, summed over the elements. With
 * log_value unset the result is the product of the F, the cdf, with no log
 * or exp of each element.
 */
template <bool upper, bool log_value = true, typename T_y, typename T_loc,
          typename T_scale, typename T_inv_scale>
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
  // to_ref only reduces compile memory by evaluating expression arguments.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));
  const bool y_inf = exp_mod_normal_check_args(function, y_val, mu_val,
                                               sigma_val, lambda_val);

  // The cdf of no elements is 1 and F is 0 or 1 at y = -inf or inf.
  constexpr double zero_value = log_value ? NEGATIVE_INFTY : 0.0;
  constexpr double one_value = log_value ? 0.0 : 1.0;
  if (size_zero(y_ref, mu_ref, sigma_ref, lambda_ref)) {
    return one_value;
  } else if (unlikely(y_inf)) {
    if (any(y_val == NEGATIVE_INFTY)) {
      return upper ? 0.0 : zero_value;
    }
    return upper ? NEGATIVE_INFTY : one_value;
  }

  scalar_seq_view<std::decay_t<decltype(y_val)>> y_vec(y_val);
  scalar_seq_view<std::decay_t<decltype(mu_val)>> mu_vec(mu_val);
  scalar_seq_view<std::decay_t<decltype(sigma_val)>> sigma_vec(sigma_val);
  scalar_seq_view<std::decay_t<decltype(lambda_val)>> lambda_vec(lambda_val);
  const std::size_t N = max_size(y_ref, mu_ref, sigma_ref, lambda_ref);
  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);
  T_partials lp = one_value;
  for (std::size_t n = 0; n < N; ++n) {
    const T_partials sigma_n = sigma_vec[n];
    const T_partials z = (y_vec[n] - mu_vec[n]) / sigma_n;
    const T_partials v = lambda_vec[n] * sigma_n;
    const auto r
        = exp_mod_normal_lcdf_kernel<upper, any_autodiff, log_value>(z, v);
    if constexpr (log_value) {
      lp += r.value;
    } else {
      lp *= r.value;
    }
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
  if constexpr (!log_value) {
    // A cdf that underflows to 0 has zero partials; the partials of log F
    // can be infinite there.
    if (unlikely(lp == 0)) {
      return make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref)
          .build(0.0);
    }
    // The partials so far are those of log F; those of the cdf scale by it.
    if constexpr (is_autodiff_v<T_y>) {
      for (std::size_t n = 0; n < stan::math::size(y_ref); ++n) {
        partials<0>(ops_partials)[n] *= lp;
      }
    }
    if constexpr (is_autodiff_v<T_loc>) {
      for (std::size_t n = 0; n < stan::math::size(mu_ref); ++n) {
        partials<1>(ops_partials)[n] *= lp;
      }
    }
    if constexpr (is_autodiff_v<T_scale>) {
      for (std::size_t n = 0; n < stan::math::size(sigma_ref); ++n) {
        partials<2>(ops_partials)[n] *= lp;
      }
    }
    if constexpr (is_autodiff_v<T_inv_scale>) {
      for (std::size_t n = 0; n < stan::math::size(lambda_ref); ++n) {
        partials<3>(ops_partials)[n] *= lp;
      }
    }
  }
  return ops_partials.build(lp);
}

}  // namespace internal

/** \ingroup prob_dists
 * Returns the log of the exponentially modified normal cumulative
 * distribution function at y, given location mu, scale sigma and inverse
 * scale lambda. Vectorized arguments give the sum over the elements.
 *
 * @tparam T_y type of the random variable
 * @tparam T_loc type of the location parameter
 * @tparam T_scale type of the scale parameter
 * @tparam T_inv_scale type of the inverse scale parameter
 * @param y random variable
 * @param mu location parameter
 * @param sigma scale parameter
 * @param lambda inverse scale parameter
 * @return log of the cumulative distribution function
 * @throw std::domain_error if y is nan, mu is infinite, or sigma or lambda is
 * not positive and finite
 * @throw std::invalid_argument if container sizes mismatch
 */
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
