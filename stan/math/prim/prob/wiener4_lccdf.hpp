#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCCDF_HPP

// #include <stan/math/prim/fun.hpp>
// #include <stan/math/prim/prob/wiener5_lpdf.hpp>
#include <stan/math/prim/prob/wiener4_lcdf.hpp>

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
 * @param w_value The relative starting point
 * @param v_value The drift rate
 * @return log probability to reach the upper bound
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob(const T_a& a, const T_v& v_value,
                        const T_w& w_value) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto v = -v_value;
  const auto w = 1 - w_value;
  if (fabs(v) == 0.0) {
    return ret_t(log1p(-w));
  }

  const auto exponent = -2.0 * v * a * (1.0 - w);
  if (exponent < 0) {
    return ret_t(log1m_exp(exponent) - log_diff_exp(2 * v * a * w, exponent));
  } else {
    return ret_t(log1m_exp(-exponent) - log1m_exp(2 * v * a));
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
 * @param w_value The relative starting point
 * @param v_value The drift rate
 * @return 'ans' term
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob_derivative_term(const T_a& a, const T_v& v_value,
                                        const T_w& w_value) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto exponent_m1 = log1p(-1.1 * 1.0e-8);
  ret_t ans;
  const auto v = -v_value;
  const auto w = 1 - w_value;

  if (fabs(v) == 0.0) {
    return ret_t(-w);
  }
  if (v < 0) {
    const auto exponent_with_1mw = 2.0 * v * a * (1.0 - w);
    const auto exponent_with_w = 2 * a * v * w;
    const auto exponent = 2 * a * v;

    if (((exponent_with_1mw >= exponent_m1) || (exponent_with_w >= exponent_m1))
        || (exponent >= exponent_m1)) {
      return ret_t(-w);
    }
    ans = LOG_TWO + exponent_with_1mw - log1m_exp(exponent_with_1mw);
    const auto diff_term = log1m_exp(exponent_with_w) - log1m_exp(exponent);
    const auto log_w = log(w);
    if (log_w > diff_term) {
      ans += log_diff_exp(log_w, diff_term);
      ans = exp(ans);
    } else {
      ans += log_diff_exp(diff_term, log_w);
      ans = -exp(ans);
    }
  } else {
    const auto exponent_with_1mw = -2.0 * v * a * (1.0 - w);
    const auto exponent = (-2 * a * v);
    if ((exponent_with_1mw >= exponent_m1) || (exponent >= exponent_m1)) {
      return ret_t(-w);
    }
    ans = LOG_TWO - log1m_exp(exponent_with_1mw);
    const auto diff_term
        = log_diff_exp(exponent_with_1mw, exponent) - log1m_exp(exponent);
    const auto log_w = log(w);
    if (log_w > diff_term) {
      ans += log_diff_exp(log_w, diff_term);
      ans = -exp(ans);
    } else {
      ans += log_diff_exp(diff_term, log_w);
      ans = exp(ans);
    }
  }
  if (fabs(ans) < INFTY) {
    return ans;
  } else {
    return ret_t(NEGATIVE_INFTY);
  }
  return ans;
}

/**
 * Calculate the derivative of the wiener probability w.r.t. 'a' (on log-scale)
 *
 * @tparam T_a type of boundary
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param a The boundary separation
 * @param w The relative starting point
 * @param v The drift rate
 * @return Gradient w.r.t. a
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob_grad_a(const T_a& a, const T_v& v,
                               const T_w& w) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  if (fabs(v) == 0.0) {
    return ret_t(0.0);
  }
  const auto deriv_term = wiener_prob_derivative_term(a, v, w);
  const auto ans = -deriv_term * v;
  if (fabs(ans) < INFTY) {
    return ans;
  } else {
    return (ret_t(NEGATIVE_INFTY));
  }
  return ans;
}

/**
 * Calculate the derivative of the wiener probability w.r.t. 'v' (on log-scale)
 *
 * @tparam T_a type of boundary
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param a The boundary separation
 * @param w The relative starting point
 * @param v The drift rate
 * @return Gradient w.r.t. v
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob_grad_v(const T_a& a, const T_v& v,
                               const T_w& w) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto deriv_term = wiener_prob_derivative_term(a, v, w);
  const auto ans = -1 * deriv_term * a;
  if (fabs(ans) < INFTY) {
    return (ans);
  } else {
    return (ret_t(NEGATIVE_INFTY));
  }
  return ans;
}

/**
 * Calculate the derivative of the wiener probability w.r.t. 'w' (on log-scale)
 *
 * @tparam T_a type of boundary
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param a The boundary separation
 * @param w_value The relative starting point
 * @param v_value The drift rate
 * @return Gradient w.r.t. w
 */
template <typename T_a, typename T_w, typename T_v>
inline auto wiener_prob_grad_w(const T_a& a, const T_v& v_value,
                               const T_w& w_value) noexcept {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  const auto v = -v_value;
  const auto w = 1 - w_value;
  if (fabs(v) == 0.0) {
    return ret_t(1 / (1.0 - w));
  }

  if (v < 0) {
    const auto exponent = 2.0 * v * a * (1.0 - w);
    return exp(LOG_TWO + exponent + log(fabs(v)) + log(a)
               - log1m_exp(exponent));
  } else {
    const auto exponent = -2.0 * v * a * (1.0 - w);
    return exp(LOG_TWO + log(fabs(v)) + log(a) - log1m_exp(exponent));
  }
}

/**
 * Calculate wiener4 ccdf (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param wildcard This parameter space is needed for a functor. Could be
 * deleted when another solution is found
 * @param err The log error tolerance
 * @return ccdf
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_wildcard, typename T_err>
inline auto wiener4_ccdf(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                         T_wildcard&& wildcard = 0.0,
                         T_err&& err = log(1e-12)) noexcept {
  const auto prob = exp(wiener_prob(a, v, w));
  const auto cdf
      = internal::wiener4_distribution<GradientCalc::ON>(y, a, v, w, 0, err);
  return prob - cdf;
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'a' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param err The log error tolerance
 * @return Gradient w.r.t. a
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_a(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err&& err = log(1e-12)) noexcept {
  const auto prob = wiener_prob(a, v, w);
  const auto prob_grad_a = wiener_prob_grad_a(a, v, w);
  const auto cdf_grad_a = wiener4_cdf_grad_a(y, a, v, w, cdf, err);
  return prob_grad_a * exp(prob) - cdf_grad_a;
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'v' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param err The log error tolerance
 * @return Gradient w.r.t. v
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_v(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err&& err = log(1e-12)) noexcept {
  const auto prob
      = wiener_prob(a, v, w);  // maybe hand over to this function, but then
                               // wiener7_integrate_cdf has problems
  const auto prob_grad_v = wiener_prob_grad_v(a, v, w);
  const auto cdf_grad_v = wiener4_cdf_grad_v(y, a, v, w, cdf, err);
  return prob_grad_v * exp(prob) - cdf_grad_v;
}

/**
 * Calculate derivative of the wiener4 ccdf w.r.t. 'w' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The CDF value
 * @param err The log error tolerance
 * @return Gradient w.r.t. w
 */
template <typename T_y, typename T_a, typename T_w, typename T_v,
          typename T_cdf, typename T_err>
inline auto wiener4_ccdf_grad_w(const T_y& y, const T_a& a, const T_v& v,
                                const T_w& w, T_cdf&& cdf,
                                T_err&& err = log(1e-12)) noexcept {
  const auto prob
      = wiener_prob(a, v, w);  // maybe hand over to this function, but then
                               // wiener7_integrate_cdf has problems
  const auto prob_grad_w = wiener_prob_grad_w(a, v, w);
  const auto cdf_grad_w = wiener4_cdf_grad_w(y, a, v, w, cdf, err);
  return prob_grad_w * exp(prob) - cdf_grad_w;
}

}  // namespace internal

/**
 * Log-CCDF for the 4-parameter Wiener distribution.
 * See 'wiener_full_lpdf' for more comprehensive documentation
 *
 * @tparam T_y type of scalar
 * @tparam T_a type of boundary
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 *
 * @param y A scalar variable; the reaction time in seconds
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
inline auto wiener_lccdf(const T_y& y, const T_a& a, const T_t0& t0,
                         const T_w& w, const T_v& v,
                         const double& precision_derivatives) {
  using T_partials_return = partials_return_t<T_y, T_a, T_t0, T_w, T_v>;
  using ret_t = return_type_t<T_y, T_a, T_t0, T_w, T_v>;

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v>::value) {
    return ret_t(0.0);
  }

  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_a_ref = ref_type_if_t<!is_constant<T_a>::value, T_a>;
  using T_t0_ref = ref_type_if_t<!is_constant<T_t0>::value, T_t0>;
  using T_w_ref = ref_type_if_t<!is_constant<T_w>::value, T_w>;
  using T_v_ref = ref_type_if_t<!is_constant<T_v>::value, T_v>;

  static constexpr const char* function_name = "wiener4_lccdf";
  if (size_zero(y, a, t0, w, v)) {
    return ret_t(0.0);
  }

  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) a_val = to_ref(as_value_column_array_or_scalar(a_ref));
  decltype(auto) v_val = to_ref(as_value_column_array_or_scalar(v_ref));
  decltype(auto) w_val = to_ref(as_value_column_array_or_scalar(w_ref));
  decltype(auto) t0_val = to_ref(as_value_column_array_or_scalar(t0_ref));
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

  const auto log_error_cdf = log(1e-6);
  const auto log_error_derivative = log(precision_derivatives);
  const T_partials_return log_error_absolute = log(1e-12);
  T_partials_return lccdf = 0.0;
  auto ops_partials
      = make_partials_propagator(y_ref, a_ref, t0_ref, w_ref, v_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;

  // calculate distribution and partials
  for (size_t i = 0; i < N; i++) {
    const auto y_value = y_vec.val(i);
    const auto a_value = a_vec.val(i);
    const auto t0_value = t0_vec.val(i);
    const auto w_value = w_vec.val(i);
    const auto v_value = v_vec.val(i);

    using internal::GradientCalc;
    const T_partials_return cdf
        = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                            GradientCalc::OFF>(
            [](auto&&... args) {
              return internal::wiener4_distribution<GradientCalc::ON>(args...);
            },
            log_error_cdf - LOG_TWO, y_value - t0_value, a_value, v_value,
            w_value, 0.0, log_error_absolute);

    const auto prob = exp(internal::wiener_prob(a_value, v_value, w_value));
    const auto ccdf = prob - cdf;

    lccdf += log(ccdf);

    const auto new_est_err = log(ccdf) + log_error_derivative - LOG_FOUR;

    const auto deriv_y = internal::estimate_with_err_check<5, 0>(
        [](auto&&... args) {
          return internal::wiener5_density<GradientCalc::ON>(args...);
        },
        new_est_err, y_value - t0_value, a_value, v_value, w_value, 0.0,
        log_error_absolute);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[i] = -deriv_y / ccdf;
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
    if (!is_constant_all<T_t0>::value) {
      partials<2>(ops_partials)[i] = deriv_y / ccdf;
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
