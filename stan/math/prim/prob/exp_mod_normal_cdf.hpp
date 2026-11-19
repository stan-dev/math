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
#include <cstddef>
#include <utility>

namespace stan {
namespace math {
namespace internal {

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


/**
 * The result from F = g Delta, Delta = E(u) - E(t) with t >= 0, for F <= 1/2:
 * F alone when calc_grad is unset, otherwise F with the partials of log F.
 * The partials of F are g times (v E(t), -E'(t) / sqrt 2,
 * v sqrt 2 (u E(t) - E'(t) / 2)), and dividing by F leaves 1 / Delta times
 * the bracket. Where 1 / Delta overflows or E'(t) underflows
 * (exp_mod_normal_use_ratio), the partials come from
 * exp_mod_normal_lcdf_terms_from_ratio.
 */
template <bool calc_grad, typename T1, typename T2, typename T3, typename T4,
          typename T5, typename T6>
inline auto exp_mod_normal_cdf_from_delta(const T1& u, const T2& t, const T3& v,
                                          const T4& log_g, const T5& e_t,
                                          const T6& dl) {
  if constexpr (!calc_grad) {
    return exp(log_g) * dl.delta;
  } else {
    const auto value = exp(log_g) * dl.delta;
    if (unlikely(exp_mod_normal_use_ratio(dl.delta, t))) {
      return exp_mod_normal_lcdf_terms_from_ratio(value, u, t, v, e_t, dl);
    }
    const auto scale = inv(dl.delta);
    const auto de_t = erfcx_derivative(t, e_t);
    return exp_mod_normal_lcdf_terms{
        value, scale * v * e_t, -scale * de_t * INV_SQRT_TWO,
        scale * v * SQRT_TWO * (u * e_t - 0.5 * de_t)};
  }
}

/**
 * F for one element: F alone when calc_grad is unset, otherwise F with the
 * partials of log F with respect to z, v and sigma d/dsigma. The regions and
 * the formulas are those of exp_mod_normal_lcdf_kernel, without the log of
 * each value.
 */
template <bool calc_grad, typename T1, typename T2>
inline auto exp_mod_normal_cdf_kernel(T1&& z, T2&& v) {
  using T_return = return_type_t<T1, T2>;
  const auto u = -z * INV_SQRT_TWO;
  const auto d = v * INV_SQRT_TWO;
  const auto t = u + d;
  const auto log_g = LOG_HALF - (0.5 * z) * z;
  if (z <= 0) {
    // F <= Phi(z) <= 1/2.
    const auto e_t = erfcx(t);
    const auto dl = exp_mod_normal_erfcx_diff<false, false, calc_grad>(
        u, t, d, erfcx(u), e_t);
    return exp_mod_normal_cdf_from_delta<calc_grad>(u, t, v, log_g, e_t, dl);
  }
  const auto log_b = d * (2 * u + d);
  // Where b_only holds, Phi(-z), Phi(v - z) and phi(z) are below 3e-19 of
  // b = exp(d (2u + d)), and b < 3/4: F = 1 - b, and the partials of F are
  // b (v, z - v, -v^2).
  const bool b_only = t < -6.5 && log_b < -0.2877;
  if constexpr (calc_grad) {
    if (b_only) {
      const auto b = exp(log_b);
      const auto r = b / (1 - b);
      return exp_mod_normal_lcdf_terms{T_return(1 - b), T_return(v * r),
                                       T_return((z - v) * r),
                                       T_return(-v * v * r)};
    }
  } else {
    if (b_only) {
      return T_return(1 - exp(log_b));
    }
  }
  const auto g = exp(log_g);
  const auto e_mu = erfcx(-u);
  const auto q = g * e_mu;
  // Assigned in one of the two branches below, so these keep a named type.
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
  const auto p = q + b;
  // F < 1/4: 1 - p would lose digits.
  const bool small_f = p > 0.75;
  if (small_f && t >= 0) {
    const auto dl = exp_mod_normal_erfcx_diff<false, false, calc_grad>(
        u, t, d, erfcx(u), e_t);
    return exp_mod_normal_cdf_from_delta<calc_grad>(u, t, v, log_g, e_t, dl);
  }
  T_return f;
  if (small_f) {
    // Here t < 0: F = -expm1(d (2u + d)) + g (E(-t) - E(-u)).
    const auto dl
        = exp_mod_normal_erfcx_diff<true, false>(-t, -u, d, e_mt, e_mu);
    f = -expm1(log_b) + g * dl.delta;
  } else {
    // F >= 1/4: F = 1 - p loses at most 2 bits.
    f = 1 - p;
  }
  if constexpr (!calc_grad) {
    return f;
  } else {
    // The partials of F: f_z, f_v, f_s.
    const T_return f_z = v * b;
    T_return f_v;
    T_return f_s;
    if (t >= 0) {
      const auto de_t = erfcx_derivative(t, e_t);
      f_v = -g * de_t * INV_SQRT_TWO;
      f_s = v * g * SQRT_TWO * (u * e_t - 0.5 * de_t);
    } else {
      const auto phi = g * SQRT_TWO_OVER_SQRT_PI;
      f_v = phi - (v - z) * b;
      f_s = v * (phi - v * b);
    }
    return exp_mod_normal_lcdf_terms{f, f_z / f, f_v / f, f_s / f};
  }
}

/**
 * The cdf: the product of the F over the elements, with no log or exp of
 * each element. The kernel returns F and the partials of log F, which are
 * scaled by the product at the end.
 */
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
  // to_ref only reduces compile memory by evaluating expression arguments.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));
  const bool y_inf = exp_mod_normal_check_args(function, y_val, mu_val,
                                               sigma_val, lambda_val);

  // The cdf of no elements is 1 and F is 0 or 1 at y = -inf or inf.
  if (size_zero(y_ref, mu_ref, sigma_ref, lambda_ref)) {
    return 1.0;
  } else if (unlikely(y_inf)) {
    return any(y_val == NEGATIVE_INFTY) ? 0.0 : 1.0;
  }

  scalar_seq_view<std::decay_t<decltype(y_val)>> y_vec(y_val);
  scalar_seq_view<std::decay_t<decltype(mu_val)>> mu_vec(mu_val);
  scalar_seq_view<std::decay_t<decltype(sigma_val)>> sigma_vec(sigma_val);
  scalar_seq_view<std::decay_t<decltype(lambda_val)>> lambda_vec(lambda_val);
  const std::size_t N = max_size(y_ref, mu_ref, sigma_ref, lambda_ref);
  T_partials cdf = 1.0;
  if constexpr (!any_autodiff) {
    // Without gradients the kernel returns the value alone.
    for (std::size_t n = 0; n < N; ++n) {
      const auto sigma_n = sigma_vec[n];
      const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
      const auto v = lambda_vec[n] * sigma_n;
      cdf *= exp_mod_normal_cdf_kernel<false>(z, v);
    }
    return cdf;
  } else {
    auto ops_partials
        = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);
    for (std::size_t n = 0; n < N; ++n) {
      const auto sigma_n = sigma_vec[n];
      const auto z = (y_vec[n] - mu_vec[n]) / sigma_n;
      const auto v = lambda_vec[n] * sigma_n;
      const auto r = exp_mod_normal_cdf_kernel<true>(z, v);
      cdf *= r.value;
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
    // A cdf that underflows to 0 has zero partials; the partials of log F
    // can be infinite there.
    if (unlikely(cdf == 0)) {
      return make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref)
          .build(0.0);
    }
    // The partials so far are those of log F; those of the cdf scale by it.
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
