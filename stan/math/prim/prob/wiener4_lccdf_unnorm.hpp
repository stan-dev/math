#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCCDF_UNNORM_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCCDF_UNNORM_HPP

#include <stan/math/prim/prob/wiener4_lcdf_unnorm.hpp>
#include <array>

namespace stan {
namespace math {
namespace internal {

/**
 * Log of probability of reaching the upper bound in diffusion process
 *
 * @tparam T_a type of boundary
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param a The boundary separation
 * @param w The relative starting point
 * @param v The drift rate
 * @return log probability to reach the upper bound
 */
template <typename T_a, typename T_w, typename T_v>
inline auto log_wiener_prob_hit_upper(const T_a& a, const T_v& v,
                                      const T_w& w) {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto neg_v = -v;
  const auto one_m_w = 1.0 - w;
  if (fabs(v) == 0.0) {
    return ret_t(log(w));
  }
  const auto exponent = 2.0 * v * a * w;
  // This branch is for numeric stability
  if (exponent < 0) {
    return ret_t(log1m_exp(exponent)
                 - log_diff_exp(2.0 * neg_v * a * one_m_w, exponent));
  } else {
    return ret_t(log1m_exp(-exponent) - log1m_exp(2.0 * neg_v * a));
  }
}

/**
 * Calculate parts of the partial derivatives for wiener_prob_grad_a and
 * wiener_prob_grad_v (on log-scale)
 *
 * @tparam T_a type of boundary
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param a The boundary separation
 * @param w The relative starting point
 * @param v The drift rate
 * @return 'ans' term
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob_derivative_term(const T_a& a, const T_v& v,
                                        const T_w& w) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto exponent_m1 = log1m(1.1 * 1.0e-8);
  const auto neg_v = -v;
  const auto one_m_w = 1 - w;
  int sign_v = neg_v < 0 ? 1 : -1;
  const auto two_a_neg_v = 2.0 * a * neg_v;
  const auto exponent_with_1mw = sign_v * two_a_neg_v * w;
  const auto exponent = sign_v * two_a_neg_v;
  const auto exponent_with_w = two_a_neg_v * one_m_w;
  // truncating longer calculations, for numerical stability
  if (unlikely((exponent_with_1mw >= exponent_m1)
               || ((exponent_with_w >= exponent_m1) && (sign_v == 1))
               || (exponent >= exponent_m1) || neg_v == 0)) {
    return ret_t(-one_m_w);
  }
  ret_t ans;
  ret_t diff_term;
  const auto log_w = log(one_m_w);
  if (neg_v < 0) {
    ans = LOG_TWO + exponent_with_1mw - log1m_exp(exponent_with_1mw);
    diff_term = log1m_exp(exponent_with_w) - log1m_exp(exponent);
  } else /* neg_v > 0 */ {
    ans = LOG_TWO - log1m_exp(exponent_with_1mw);
    diff_term = log_diff_exp(exponent_with_1mw, exponent) - log1m_exp(exponent);
  }
  if (log_w > diff_term) {
    ans = sign_v * exp(ans + log_diff_exp(log_w, diff_term));
  } else {
    ans = -sign_v * exp(ans + log_diff_exp(diff_term, log_w));
  }
  if (unlikely(!is_scal_finite(ans))) {
    return ret_t(NEGATIVE_INFTY);
  }
  return ans;
}

/**
 * Log of the wiener4 ccdf, P - F, from its large-time series, and the
 * derivatives of this log with respect to a, w and v. With u = 1 - w and
 * b = -v (upper boundary),
 *   P - F = 2 pi / a^2 exp(-u a b - b^2 y / 2) S,
 *   S = sum_{k >= 1} k sin(k pi u) e_k / d_k,
 *   e_k = exp(-k^2 pi^2 y / (2 a^2)),  d_k = b^2 + k^2 pi^2 / a^2.
 * Where F is close to P, the difference P - F loses digits; the series
 * does not. The sum stops when a bound of the remaining terms (and, if
 * Grad, of the remaining terms of the derivative sums) is below
 * exp(log_err) |S|. The bound uses |k sin(k pi u) e_k / d_k| <= a^2 / pi^2
 * e_k / k.
 *
 * @tparam Grad Whether the derivatives are computed
 * @tparam T_y type of reaction time
 * @tparam T_a type of boundary separation
 * @tparam T_v type of drift rate
 * @tparam T_w type of relative starting point
 * @tparam T_err type of log error tolerance
 *
 * @param y The reaction time in seconds, greater than 0
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param log_err The log error tolerance, relative to P - F
 * @return log(P - F), then its derivatives with respect to a, w and v
 * (0 if not Grad)
 */
template <bool Grad, typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_err>
inline auto wiener4_log_ccdf_series(const T_y& y, const T_a& a, const T_v& v,
                                    const T_w& w, T_err log_err) {
  using ret_t = return_type_t<T_y, T_a, T_v, T_w>;
  const auto u = 1.0 - w;
  const auto b = -v;
  const auto pi_a_sq = square(pi() / a);
  const auto c = 0.5 * pi_a_sq * y;
  const auto log_a_pi_sq = 2.0 * (log(a) - LOG_PI);
  // bound of |derivative term| / (k |value term bound|)
  const auto log_grad_mult
      = log(pi_a_sq * y / a + 2.0 / a + pi() + 2.0 * fabs(b) / pi_a_sq);
  ret_t s_pos = NEGATIVE_INFTY;
  ret_t s_neg = NEGATIVE_INFTY;
  ret_t ga_pos = NEGATIVE_INFTY;
  ret_t ga_neg = NEGATIVE_INFTY;
  ret_t gu_pos = NEGATIVE_INFTY;
  ret_t gu_neg = NEGATIVE_INFTY;
  ret_t gb_pos = NEGATIVE_INFTY;
  ret_t gb_neg = NEGATIVE_INFTY;
  const auto add
      = [](ret_t& pos, ret_t& neg, const ret_t& log_term, bool positive) {
          if (positive) {
            pos = log_sum_exp(pos, log_term);
          } else {
            neg = log_sum_exp(neg, log_term);
          }
        };
  constexpr int max_terms = 100000;
  for (int k = 1; k <= max_terms; ++k) {
    const double kd = k;
    const auto d_k = square(b) + kd * kd * pi_a_sq;
    const auto log_e_d = -kd * kd * c - log(d_k);
    const auto sin_k = sin(kd * pi() * u);
    if (sin_k != 0) {
      const ret_t log_term = log(kd) + log_e_d + log(fabs(sin_k));
      add(s_pos, s_neg, log_term, sin_k > 0);
      if constexpr (Grad) {
        add(ga_pos, ga_neg,
            log_term + log(kd * kd * pi_a_sq * (y + 2.0 / d_k) / a), sin_k > 0);
        if (b != 0) {
          add(gb_pos, gb_neg, log_term + log(2.0 * fabs(b) / d_k),
              (sin_k > 0) == (b < 0));
        }
      }
    }
    if constexpr (Grad) {
      const auto cos_k = cos(kd * pi() * u);
      if (cos_k != 0) {
        add(gu_pos, gu_neg, 2.0 * log(kd) + LOG_PI + log_e_d + log(fabs(cos_k)),
            cos_k > 0);
      }
    }
    const auto ratio = (Grad ? 2.0 : 1.0) * exp(-(2.0 * kd + 3.0) * c);
    if (s_pos > s_neg && ratio < 0.5) {
      const auto log_tail
          = log_a_pi_sq - square(kd + 1.0) * c - log1m(ratio)
            + (Grad ? log(kd + 1.0) + log_grad_mult : -log(kd + 1.0));
      if (log_tail <= log_err + log_diff_exp(s_pos, s_neg)) {
        break;
      }
    }
  }
  const ret_t log_s = s_pos > s_neg ? ret_t(log_diff_exp(s_pos, s_neg))
                                    : ret_t(NEGATIVE_INFTY);
  std::array<ret_t, 4> out{ret_t(LOG_TWO + LOG_PI - 2.0 * log(a) - u * a * b
                                 - 0.5 * square(b) * y + log_s),
                           ret_t(0.0), ret_t(0.0), ret_t(0.0)};
  if constexpr (Grad) {
    const auto ratio_to_s = [&log_s](const ret_t& pos, const ret_t& neg) {
      if (pos > neg) {
        return ret_t(exp(log_diff_exp(pos, neg) - log_s));
      }
      if (neg > pos) {
        return ret_t(-exp(log_diff_exp(neg, pos) - log_s));
      }
      return ret_t(0.0);
    };
    out[1] = -2.0 / a - u * b + ratio_to_s(ga_pos, ga_neg);
    out[2] = a * b - ratio_to_s(gu_pos, gu_neg);
    out[3] = u * a + b * y - ratio_to_s(gb_pos, gb_neg);
  }
  return out;
}

/**
 * Log of the wiener4 ccdf, P - F: log(P - F) where F < P / 2, otherwise the
 * large-time series of P - F (wiener4_log_ccdf_series). For y <= 0 the
 * ccdf is P.
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param log_err The log error tolerance, relative to P (for F) and to
 * P - F (for the series)
 * @return log ccdf
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_err>
inline auto wiener4_log_ccdf(const T_y& y, const T_a& a, const T_v& v,
                             const T_w& w, T_err log_err) {
  using ret_t = return_type_t<T_y, T_a, T_v, T_w>;
  const auto log_p = log_wiener_prob_hit_upper(a, v, w);
  if (y <= 0) {
    return ret_t(log_p);
  }
  const auto log_f = wiener4_distribution<GradientCalc::OFF>(
      y, a, v, w, log_err + fmin(0.0, fmax(log_p, -1e3)));
  if (log_f < log_p - LOG_TWO) {
    return ret_t(log_diff_exp(log_p, log_f));
  }
  return ret_t(wiener4_log_ccdf_series<false>(y, a, v, w, log_err)[0]);
}

/**
 * Calculate wiener4 ccdf (natural-scale), see wiener4_log_ccdf
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return ccdf
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_err>
inline auto wiener4_ccdf(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                         T_err log_err = log(1e-12)) noexcept {
  return exp(wiener4_log_ccdf(y, a, v, w, log_err));
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'a' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to a
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_a(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err log_err = log(1e-12)) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;

  using res_t = return_type_t<T_y, T_a, T_w, T_v, T_cdf>;
  // derivative of the wiener probability w.r.t. 'a' (on log-scale)
  auto prob_grad_a = -wiener_prob_derivative_term(a, v, w) * v;
  if (!is_scal_finite(prob_grad_a)) {
    prob_grad_a = ret_t(NEGATIVE_INFTY);
  }
  const auto log_prob_hit_upper = log_wiener_prob_hit_upper(a, v, w);
  if (y <= 0) {
    return res_t(prob_grad_a * exp(log_prob_hit_upper));
  }
  if (log(cdf) >= log_prob_hit_upper - LOG_TWO) {
    const auto series = wiener4_log_ccdf_series<true>(y, a, v, w, log_err);
    return res_t(exp(series[0]) * series[1]);
  }
  const auto cdf_grad_a = wiener4_cdf_grad_a(y, a, v, w, cdf, log_err);
  return res_t(prob_grad_a * exp(log_prob_hit_upper) - cdf_grad_a);
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'v' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to v
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_v(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err log_err = log(1e-12)) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  using res_t = return_type_t<T_y, T_a, T_w, T_v, T_cdf>;
  const auto log_prob_hit_upper = log_wiener_prob_hit_upper(a, v, w);
  // derivative of the wiener probability w.r.t. 'v' (on log-scale)
  auto prob_grad_v = -wiener_prob_derivative_term(a, v, w) * a;
  if (!is_scal_finite(fabs(prob_grad_v))) {
    prob_grad_v = ret_t(NEGATIVE_INFTY);
  }
  if (y <= 0) {
    return res_t(prob_grad_v * exp(log_prob_hit_upper));
  }
  if (log(cdf) >= log_prob_hit_upper - LOG_TWO) {
    const auto series = wiener4_log_ccdf_series<true>(y, a, v, w, log_err);
    return res_t(exp(series[0]) * series[3]);
  }
  const auto cdf_grad_v = wiener4_cdf_grad_v(y, a, v, w, cdf, log_err);
  return res_t(prob_grad_v * exp(log_prob_hit_upper) - cdf_grad_v);
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'w' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to w
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_w(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err log_err = log(1e-12)) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto log_prob_hit_upper = log_wiener_prob_hit_upper(a, v, w);
  // derivative of the wiener probability w.r.t. 'v' (on log-scale)
  const auto exponent = -sign(v) * 2.0 * v * a * w;
  auto prob_grad_w
      = (v != 0) ? exp(LOG_TWO + log(fabs(v)) + log(a) - log1m_exp(exponent))
                 : ret_t(1 / w);
  if (v > 0) {
    prob_grad_w *= exp(exponent);
  }
  using res_t = return_type_t<T_y, T_a, T_w, T_v, T_cdf>;
  if (y <= 0) {
    return res_t(prob_grad_w * exp(log_prob_hit_upper));
  }
  if (log(cdf) >= log_prob_hit_upper - LOG_TWO) {
    const auto series = wiener4_log_ccdf_series<true>(y, a, v, w, log_err);
    return res_t(exp(series[0]) * series[2]);
  }
  const auto cdf_grad_w = wiener4_cdf_grad_w(y, a, v, w, cdf, log_err);
  return res_t(prob_grad_w * exp(log_prob_hit_upper) - cdf_grad_w);
}

}  // namespace internal

/**
 * Log-CCDF for the 4-parameter Wiener distribution.
 * See 'wiener_full_lpdf' for more comprehensive documentation.
 *
 * @tparam T_y type of reaction time
 * @tparam T_a type of boundary
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param t0 The non-decision time
 * @param w The relative starting point
 * @param v The drift rate
 * @param precision_derivatives Level of precision in estimation
 * @return The log of the Wiener first passage time distribution with
 *  the specified arguments for upper boundary responses
 */
template <bool propto = false, typename T_y, typename T_a, typename T_t0,
          typename T_w, typename T_v>
inline auto wiener_lccdf_unnorm(const T_y& y, const T_a& a, const T_t0& t0,
                                const T_w& w, const T_v& v,
                                const double& precision_derivatives = 1e-4) {
  using T_partials_return = partials_return_t<T_y, T_a, T_t0, T_w, T_v>;
  using ret_t = return_type_t<T_y, T_a, T_t0, T_w, T_v>;
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_v_ref = ref_type_t<T_v>;
  using internal::GradientCalc;

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;

  auto y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  auto a_val = to_ref(as_value_column_array_or_scalar(a_ref));
  auto v_val = to_ref(as_value_column_array_or_scalar(v_ref));
  auto w_val = to_ref(as_value_column_array_or_scalar(w_ref));
  auto t0_val = to_ref(as_value_column_array_or_scalar(t0_ref));

  static constexpr const char* function_name = "wiener4_lccdf";
  if (size_zero(y, a, t0, w, v)) {
    return ret_t(0.0);
  }

  if constexpr (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v>::value) {
    return ret_t(0.0);
  }

  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0);
  check_positive_finite(function_name, "Random variable", y_val);
  check_positive_finite(function_name, "Boundary separation", a_val);
  check_finite(function_name, "Drift rate", v_val);
  check_less(function_name, "A-priori bias", w_val, 1);
  check_greater(function_name, "A-priori bias", w_val, 0);
  check_nonnegative(function_name, "Nondecision time", t0_val);
  check_finite(function_name, "Nondecision time", t0_val);

  const size_t N = max_size(y, a, t0, w, v);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  const size_t N_y_t0 = max_size(y, t0);

  for (size_t i = 0; i < N_y_t0; ++i) {
    if (y_vec[i] <= t0_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << t0_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }

  // for precs. 1e-6, 1e-12, see Hartmann et al. (2021), Henrich et al. (2023)
  const auto log_error_derivative = log(precision_derivatives);
  const T_partials_return log_error_absolute = log(1e-12);
  T_partials_return lccdf = 0.0;
  auto ops_partials
      = make_partials_propagator(y_ref, a_ref, t0_ref, w_ref, v_ref);

  const double LOG_FOUR = std::log(4.0);

  // calculate distribution and partials
  for (size_t i = 0; i < N; i++) {
    const auto y_value = y_vec.val(i);
    const auto a_value = a_vec.val(i);
    const auto t0_value = t0_vec.val(i);
    const auto w_value = w_vec.val(i);
    const auto v_value = v_vec.val(i);

    const T_partials_return y_diff = y_value - t0_value;
    const T_partials_return log_prob_hit_upper
        = internal::log_wiener_prob_hit_upper(a_value, v_value, w_value);
    // F with a tolerance relative to P
    const T_partials_return log_cdf
        = internal::wiener4_distribution<GradientCalc::OFF>(
            y_diff, a_value, v_value, w_value,
            log_error_absolute + fmin(0.0, fmax(log_prob_hit_upper, -1e3)));
    // Where F is close to P, P - F loses digits: use the large-time series
    // of P - F there, and its derivatives.
    const bool use_series = !(log_cdf < log_prob_hit_upper - LOG_TWO);
    constexpr bool grad_series = !is_constant_all<T_a, T_w, T_v>::value;
    std::array<T_partials_return, 4> series{0.0, 0.0, 0.0, 0.0};
    T_partials_return log_ccdf_single_value;
    if (use_series) {
      series = internal::wiener4_log_ccdf_series<grad_series>(
          y_diff, a_value, v_value, w_value, log_error_absolute);
      log_ccdf_single_value = series[0];
    } else {
      log_ccdf_single_value = log_diff_exp(log_prob_hit_upper, log_cdf);
    }
    const T_partials_return cdf = exp(log_cdf);
    const T_partials_return ccdf = exp(log_ccdf_single_value);

    lccdf += log_ccdf_single_value;

    const auto new_est_err
        = log_ccdf_single_value + log_error_derivative - LOG_FOUR;

    if constexpr (!is_constant_all<T_y>::value
                  || !is_constant_all<T_t0>::value) {
      // density / ccdf, on the log scale, with a density tolerance relative
      // to the ccdf
      const auto deriv_y
          = exp(internal::wiener5_density<GradientCalc::OFF>(
                    y_diff, a_value, v_value, w_value, 0.0,
                    log_error_absolute + fmin(0.0, log_ccdf_single_value))
                - log_ccdf_single_value);
      if constexpr (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[i] = -deriv_y;
      }
      if constexpr (!is_constant_all<T_t0>::value) {
        partials<2>(ops_partials)[i] = deriv_y;
      }
    }
    if constexpr (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = use_series ? series[1]
                       : internal::estimate_with_err_check<5, 0>(
                             [](auto&&... args) {
                               return internal::wiener4_ccdf_grad_a(args...);
                             },
                             new_est_err, y_diff, a_value, v_value, w_value,
                             cdf, log_error_absolute)
                             / ccdf;
    }
    if constexpr (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = use_series ? series[2]
                       : internal::estimate_with_err_check<5, 0>(
                             [](auto&&... args) {
                               return internal::wiener4_ccdf_grad_w(args...);
                             },
                             new_est_err, y_diff, a_value, v_value, w_value,
                             cdf, log_error_absolute)
                             / ccdf;
    }
    if constexpr (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = use_series ? series[3]
                       : internal::wiener4_ccdf_grad_v(y_diff, a_value, v_value,
                                                       w_value, cdf,
                                                       log_error_absolute)
                             / ccdf;
    }
  }  // for loop
  return ops_partials.build(lccdf);
}
}  // namespace math
}  // namespace stan
#endif
