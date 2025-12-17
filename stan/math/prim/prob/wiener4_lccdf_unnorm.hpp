#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCCDF_UNNORM_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCCDF_UNNORM_HPP

#include <stan/math/prim/prob/wiener4_lcdf_unnorm.hpp>

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
 * Calculate wiener4 ccdf (natural-scale)
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
  const auto prob_hit_upper = exp(log_wiener_prob_hit_upper(a, v, w));
  const auto cdf
      = internal::wiener4_distribution<GradientCalc::ON>(y, a, v, w, log_err);
  return prob_hit_upper - cdf;
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

  // derivative of the wiener probability w.r.t. 'a' (on log-scale)
  auto prob_grad_a = -wiener_prob_derivative_term(a, v, w) * v;
  if (!is_scal_finite(prob_grad_a)) {
    prob_grad_a = ret_t(NEGATIVE_INFTY);
  }
  const auto log_prob_hit_upper = log_wiener_prob_hit_upper(a, v, w);
  const auto cdf_grad_a = wiener4_cdf_grad_a(y, a, v, w, cdf, log_err);
  return prob_grad_a * exp(log_prob_hit_upper) - cdf_grad_a;
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
  const auto log_prob_hit_upper = log_wiener_prob_hit_upper(a, v, w);
  // derivative of the wiener probability w.r.t. 'v' (on log-scale)
  auto prob_grad_v = -wiener_prob_derivative_term(a, v, w) * a;
  if (!is_scal_finite(fabs(prob_grad_v))) {
    prob_grad_v = ret_t(NEGATIVE_INFTY);
  }

  const auto cdf_grad_v = wiener4_cdf_grad_v(y, a, v, w, cdf, log_err);
  return prob_grad_v * exp(log_prob_hit_upper) - cdf_grad_v;
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

  const auto cdf_grad_w = wiener4_cdf_grad_w(y, a, v, w, cdf, log_err);
  return prob_grad_w * exp(log_prob_hit_upper) - cdf_grad_w;
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

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v>::value) {
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
  const auto log_error_cdf = log(1e-6);
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

    const T_partials_return cdf
        = internal::estimate_with_err_check<4, 0, GradientCalc::OFF,
                                            GradientCalc::OFF>(
            [](auto&&... args) {
              return internal::wiener4_distribution<GradientCalc::ON>(args...);
            },
            log_error_cdf - LOG_TWO, y_value - t0_value, a_value, v_value,
            w_value, log_error_absolute);

    const auto prob_hit_upper
        = exp(internal::log_wiener_prob_hit_upper(a_value, v_value, w_value));
    const auto ccdf = prob_hit_upper - cdf;
    const auto log_ccdf_single_value = log(ccdf);

    lccdf += log_ccdf_single_value;

    const auto new_est_err
        = log_ccdf_single_value + log_error_derivative - LOG_FOUR;

    if (!is_constant_all<T_y>::value || !is_constant_all<T_t0>::value) {
      const auto deriv_y = internal::estimate_with_err_check<5, 0>(
          [](auto&&... args) {
            return internal::wiener5_density<GradientCalc::ON>(args...);
          },
          new_est_err, y_value - t0_value, a_value, v_value, w_value, 0.0,
          log_error_absolute);
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[i] = -deriv_y / ccdf;
      }
      if (!is_constant_all<T_t0>::value) {
        partials<2>(ops_partials)[i] = deriv_y / ccdf;
      }
    }
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0>(
                [](auto&&... args) {
                  return internal::wiener4_ccdf_grad_a(args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / ccdf;
    }
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0>(
                [](auto&&... args) {
                  return internal::wiener4_ccdf_grad_w(args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / ccdf;
    }
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener4_ccdf_grad_v(y_value - t0_value, a_value, v_value,
                                          w_value, cdf, log_error_absolute)
            / ccdf;
    }
  }  // for loop
  return ops_partials.build(lccdf);
}
}  // namespace math
}  // namespace stan
#endif
