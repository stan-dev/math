#ifndef STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Calculate the 'error_term' term for a wiener5 density or gradient
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @return 'error_term' term
 */
inline double wiener5_compute_error_term(double y, double a, double v_value,
                                         double w_value, double sv) noexcept {
  const double w = 1.0 - w_value;
  const double v = -v_value;
  const double sv_sqr = square(sv);
  const double one_plus_svsqr_y = 1 + sv_sqr * y;
  const double two_avw = 2 * a * v * w;
  const double two_log_a = 2 * log(a);
  if (sv != 0) {
    return (sv_sqr * square(a * w) - two_avw - square(v) * y) / 2.0
               / one_plus_svsqr_y
           - two_log_a - 0.5 * log(one_plus_svsqr_y);
  } else {
    return (-two_avw - square(v) * y) / 2.0 - two_log_a;
  }
}

/**
 * Calculate the 'density_part_one' term for a wiener5 density or gradient
 *
 * @tparam GradA Whether the calculation is for gradient w.r.t. 'a'
 * @tparam GradT Whether the calculation is for gradient w.r.t. 't'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @return 'density_part_one' term
 */
template <bool GradA, bool GradT>
inline double wiener5_density_part_one(double y, double a, double v_value,
                                       double w_value, double sv) noexcept {
  const double w = 1.0 - w_value;
  const double v = -v_value;
  const double sv_sqr = square(sv);
  const double one_plus_svsqr_y = 1 + sv_sqr * y;
  const double two_avw = 2 * a * v * w;

  const double var_a = GradA ? w : a;
  const double var_b = GradA ? a : w;

  if (GradT) {
    if (sv != 0) {
      return -0.5
             * (square(sv_sqr) * (y + square(a * w)) + sv_sqr * (1 - two_avw)
                + square(v))
             / square(one_plus_svsqr_y);
    } else {
      return -0.5 * square(v);
    }
  }

  if (sv != 0) {
    return (-v * var_a + sv_sqr * square(var_a) * var_b) / one_plus_svsqr_y;
  } else {
    return -v * var_a;
  }
}

/**
 * Calculate the 'n_terms_small_t' term for a wiener5 density or gradient
 *
 * @tparam Density Whether the calculation is for the density
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param w_value The relative starting point
 * @param error The error tolerance
 * @return 'n_terms_small_t' term
 */
template <bool Density, bool GradW>
inline double wiener5_n_terms_small_t(double y, double a, double w_value,
                                      double error) noexcept {
  const double two_error = 2.0 * error;
  const double y_asq = y / square(a);
  const double two_log_a = 2 * log(a);
  const double log_y_asq = log(y) - two_log_a;
  const double w = 1.0 - w_value;

  const double n_1_factor = Density ? 2 : 3;
  const double n_1 = (sqrt(n_1_factor * y_asq) + w) / 2.0;
  double u_eps;
  if (Density || GradW) {
    u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log_y_asq + two_error);
  } else {
    u_eps = fmin(-3.0,
                 (log(8.0) - log(27.0) + LOG_PI + 4.0 * log_y_asq + two_error));
  }
  const double arg_mult = (Density || GradW) ? 1 : 3;
  const double arg = -arg_mult * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));

  const double n_2
      = (arg > 0) ? GradW ? 0.5 * (sqrt(arg) + w) : 0.5 * (sqrt(arg) - w) : n_1;
  return ceil(fmax(n_1, n_2));
}

/**
 * Calculate the 'n_terms_small_t' term for a wiener5 density or gradient
 *
 * @tparam Density Whether the calculation is for the density
 * @tparam GradT Whether the calculation is for gradient w.r.t. 't'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param w_value The relative starting point
 * @param error The error tolerance
 * @return 'n_terms_large_t' term
 */
template <bool Density, bool GradW>
inline double wiener5_n_terms_largel_t(double y, double a, double w_value,
                                       double error) noexcept {
  const double two_error = 2.0 * error;
  const double y_asq = y / square(a);
  const double two_log_a = 2 * log(a);
  const double log_y_asq = log(y) - two_log_a;
  const double w = 1.0 - w_value;

  const double n_1_factor = GradW ? 2 : 3;
  static constexpr double PI_SQUARED = pi() * pi();

  double n_1;
  double n_2;
  if (Density) {
    n_1 = 1.0 / (pi() * sqrt(y_asq));
    const double two_log_piy = -2.0 * (LOG_PI + log_y_asq + error);
    n_2 = (two_log_piy >= 0) ? sqrt(two_log_piy / (PI_SQUARED * y_asq)) : 0.0;
  } else {
    n_1 = sqrt(n_1_factor / y_asq) / pi();
    const double u_eps_arg
        = GradW
              ? log(4.0) - log(9.0) + 2.0 * LOG_PI + 3.0 * log_y_asq + two_error
              : log(3.0) - log(5.0) + LOG_PI + 2.0 * log_y_asq + error;
    const double u_eps = fmin(-1, u_eps_arg);
    const double arg_mult = GradW ? 1 : (2.0 / PI_SQUARED / y_asq);
    const double arg = -arg_mult * (u_eps - sqrt(-2.0 * u_eps - 2.0));
    n_2 = GradW ? (arg > 0) ? sqrt(arg / y_asq) / pi() : n_1
                : (arg > 0) ? sqrt(arg) : n_1;
  }

  return ceil(fmax(n_1, n_2));
}

/**
 * Calculate the 'result' term and its sign for a wiener5 density or gradient
 *
 * @tparam Density Whether the calculation is for the density
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param w_value The relative starting point
 * @param n_terms_small_t The n_terms_small_t term
 * @param n_terms_large_t The n_terms_large_t term
 * @return 'result' sum and its sign
 */
template <bool Density, bool GradW>
inline std::tuple<double, int> wiener5_log_sum_exp(
    double y, double a, double w_value, size_t n_terms_small_t,
    size_t n_terms_large_t) noexcept {
  const double y_asq = y / square(a);
  const double w = 1.0 - w_value;
  const bool small_n_terms_small_t
      = Density ? (2 * n_terms_small_t <= n_terms_large_t)
                : (2 * n_terms_small_t < n_terms_large_t);
  const double scaling = small_n_terms_small_t ? inv(2.0 * y_asq) : y_asq / 2.0;

  double prev_val = NEGATIVE_INFTY;
  double current_val = NEGATIVE_INFTY;
  int prev_sign = 1;
  int current_sign = 1;

  if (small_n_terms_small_t) {
    const double mult = Density ? 1 : 3;
    const double offset = GradW ? y_asq : 0;
    const double sqrt_offset = sqrt(offset);
    for (size_t k = n_terms_small_t; k >= 1; k--) {
      const double wp2k = w + 2.0 * k;
      const double wm2k = w - 2.0 * k;
      int wp2k_sign = (wp2k > sqrt_offset) ?: 1 - 1;
      int wm2k_sign = (wm2k > sqrt_offset) ? 1 : -1;
      double wp2k_quant
          = GradW ? log(fabs((square(wp2k) - offset))) - square(wp2k) * scaling
                  : mult * log(wp2k_sign * wp2k) - square(wp2k) * scaling;
      double wm2k_quant
          = GradW ? log(fabs((square(wm2k) - offset))) - square(wm2k) * scaling
                  : mult * log(wm2k_sign * wm2k) - square(wm2k) * scaling;
      double k_term;
      int k_sign;
      std::forward_as_tuple(k_term, k_sign) = log_sum_exp_signed(
          wm2k_quant, -1 * wm2k_sign, wp2k_quant, wp2k_sign);
      std::forward_as_tuple(current_val, current_sign)
          = log_sum_exp_signed(k_term, k_sign, prev_val, prev_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
    double new_val = GradW ? log(fabs(square(w) - offset)) - square(w) * scaling
                           : mult * log(w) - square(w) * scaling;
    int new_val_sign
        = GradW ? (w > sqrt_offset ? 1 : -1) : (new_val > 0 ? 1 : -1);
    int factor_sign = GradW ? 1 : -1;
    std::forward_as_tuple(current_val, current_sign)
        = log_sum_exp_signed(new_val, factor_sign * new_val_sign, current_val,
                             factor_sign * current_sign);
  } else {
    double mult = 3;
    if (Density) {
      mult = 1;
    } else if (GradW) {
      mult = 2;
    }
    for (size_t k = n_terms_large_t; k >= 1; k--) {
      const double pi_k = k * pi();
      const double check = (GradW) ? cos(pi_k * w) : sin(pi_k * w);
      int check_sign = sign(check);
      double n_terms_large_t_quant
          = mult * log(k) - square(pi_k) * scaling + log(fabs(check));
      std::forward_as_tuple(current_val, current_sign) = log_sum_exp_signed(
          prev_val, prev_sign, n_terms_large_t_quant, check_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
  }
  return std::make_tuple(current_val, current_sign);
}

/**
 * Calculate the wiener5 density
 *
 * @tparam NaturalScale Whether to return the density on natural (true) or
 * log-scale (false)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return density
 */
template <bool NaturalScale = false>
inline double wiener5_density(double y, double a, double v_value,
                              double w_value, double sv,
                              double err = log(1e-12)) noexcept {
  const double error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const double error = (err - error_term);
  const double n_terms_small_t
      = wiener5_n_terms_small_t<true, false>(y, a, w_value, error);
  const double n_terms_large_t
      = wiener5_n_terms_largel_t<true, false>(y, a, w_value, error);

  double result;
  int newsign;
  double log_density;
  std::forward_as_tuple(result, newsign) = wiener5_log_sum_exp<true, false>(
      y, a, w_value, n_terms_small_t, n_terms_large_t);
  if (2 * n_terms_small_t <= n_terms_large_t) {
    log_density = error_term - 0.5 * LOG_TWO - LOG_SQRT_PI
                  - 1.5 * (log(y) - 2 * log(a)) + result;
  } else {
    log_density = error_term + result + LOG_PI;
  }

  return NaturalScale ? exp(log_density) : log_density;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 't'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. t
 */
template <bool WrtLog = false>
inline double grad_wiener5_t(double y, double a, double v_value, double w_value,
                             double sv, double err = log(1e-12)) noexcept {
  const double two_log_a = 2 * log(a);
  const double log_y_asq = log(y) - two_log_a;
  const double error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const double density_part_one
      = wiener5_density_part_one<false, true>(y, a, v_value, w_value, sv);
  const double error = (err - error_term) + two_log_a;

  const double n_terms_small_t
      = wiener5_n_terms_small_t<false, false>(y, a, w_value, error);
  const double n_terms_large_t
      = wiener5_n_terms_largel_t<false, false>(y, a, w_value, error);
  double result;
  int newsign;
  std::forward_as_tuple(result, newsign) = wiener5_log_sum_exp<false, false>(
      y, a, w_value, n_terms_small_t, n_terms_large_t);

  const double error_log_density
      = log(max(fabs(density_part_one - 1.5 / y), fabs(density_part_one)));
  const double log_density
      = wiener5_density(y, a, v_value, w_value, sv, err - error_log_density);
  double ans;
  if (2 * n_terms_small_t < n_terms_large_t) {
    ans = density_part_one - 1.5 / y
          + newsign
                * exp(error_term - two_log_a - 1.5 * LOG_TWO - LOG_SQRT_PI
                      - 3.5 * log_y_asq + result - log_density);
  } else {
    ans = density_part_one
          - newsign
                * exp(error_term - two_log_a + 3.0 * LOG_PI - LOG_TWO + result
                      - log_density);
  }
  return WrtLog ? ans * exp(log_density) : ans;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'a'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. a
 */
template <bool WrtLog = false>
inline double grad_wiener5_a(double y, double a, double v_value, double w_value,
                             double sv, double err = log(1e-12)) noexcept {
  const double two_log_a = 2 * log(a);
  const double log_y_asq = log(y) - two_log_a;
  const double error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const double density_part_one
      = wiener5_density_part_one<true, false>(y, a, v_value, w_value, sv);
  const double error = err - error_term + 3 * log(a) - log(y) - LOG_TWO;

  const double n_terms_small_t
      = wiener5_n_terms_small_t<false, false>(y, a, w_value, error);
  const double n_terms_large_t
      = wiener5_n_terms_largel_t<false, false>(y, a, w_value, error);
  double result;
  int newsign;
  std::forward_as_tuple(result, newsign) = wiener5_log_sum_exp<false, false>(
      y, a, w_value, n_terms_small_t, n_terms_large_t);

  const double error_log_density = log(
      max(fabs(density_part_one + 1.0 / a), fabs(density_part_one - 2.0 / a)));
  const double log_density
      = wiener5_density(y, a, v_value, w_value, sv, err - error_log_density);
  double ans;
  if (2 * n_terms_small_t < n_terms_large_t) {
    ans = density_part_one + 1.0 / a
          - newsign
                * exp(-0.5 * LOG_TWO - LOG_SQRT_PI - 2.5 * log(y)
                      + 2.0 * two_log_a + error_term + result - log_density);
  } else {
    ans = density_part_one - 2.0 / a
          + newsign
                * exp(log(y) + error_term - 3 * (log(a) - LOG_PI) + result
                      - log_density);
  }
  return WrtLog ? ans * exp(log_density) : ans;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'v'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. v
 */
template <bool WrtLog = false>
inline double grad_wiener5_v(double y, double a, double v_value, double w_value,
                             double sv, double err = log(1e-12)) noexcept {
  double ans = (a * (1 - w_value) - v_value * y);
  if (sv != 0) {
    ans /= 1 + square(sv) * y;
  }
  return WrtLog ? ans * wiener5_density<true>(y, a, v_value, w_value, sv, err)
                : ans;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'w'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. w
 */
template <bool WrtLog = false>
inline double grad_wiener5_w(double y, double a, double v_value, double w_value,
                             double sv, double err = log(1e-12)) noexcept {
  const double two_log_a = 2 * log(a);
  const double log_y_asq = log(y) - two_log_a;
  const double error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const double density_part_one
      = wiener5_density_part_one<false, false>(y, a, v_value, w_value, sv);
  const double error = (err - error_term);

  const double n_terms_small_t
      = wiener5_n_terms_small_t<false, true>(y, a, w_value, error);
  const double n_terms_large_t
      = wiener5_n_terms_largel_t<false, true>(y, a, w_value, error);
  double result;
  int newsign;
  std::forward_as_tuple(result, newsign) = wiener5_log_sum_exp<false, true>(
      y, a, w_value, n_terms_small_t, n_terms_large_t);

  const double log_density = wiener5_density(y, a, v_value, w_value, sv,
                                             err - log(fabs(density_part_one)));
  double ans;
  if (2 * n_terms_small_t < n_terms_large_t) {
    ans = -(density_part_one
            - newsign
                  * exp(result - (log_density - error_term) - 2.5 * log_y_asq
                        - 0.5 * LOG_TWO - 0.5 * LOG_PI));
  } else {
    ans = -(density_part_one
            + newsign * exp(result - (log_density - error_term) + 2 * LOG_PI));
  }
  return WrtLog ? ans * exp(log_density) : ans;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'sv'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. sv
 */
template <bool WrtLog = false>
inline double grad_wiener5_sv(double y, double a, double v_value,
                              double w_value, double sv,
                              double err = log(1e-12)) noexcept {
  const double one_plus_svsqr_y = 1 + square(sv) * y;
  const double w = 1.0 - w_value;
  const double v = -v_value;
  const double t1 = -y / one_plus_svsqr_y;
  const double t2 = (square(a * w) + 2 * a * v * w * y + square(v * y))
                    / square(one_plus_svsqr_y);
  const double ans = sv * (t1 + t2);
  return WrtLog ? ans * wiener5_density<true>(y, a, v_value, w_value, sv, err)
                : ans;
}

/**
 * Utility function for replacing a value with a specified error value
 */
template <size_t NestedIndex>
inline void assign_err(double arg, double err) {
  arg = err;
}

/**
 * Utility function for replacing a value with a specified error value,
 * overload for use when the value is stored within a tuple.
 */
template <size_t NestedIndex, typename... TArgs>
inline void assign_err(std::tuple<TArgs...>& args_tuple, double err) {
  std::get<NestedIndex>(args_tuple) = err;
}

/**
 * Utility function for estimating a function with a given set of arguments,
 * checking the result against a provided error tolerance, and re-estimating
 * the function with the increased error if it fails.
 *
 * @tparam ErrIndex Position of the error argument in the provided arguments
 * @tparam NestedIndex Nested position if the error argument is within a tuple
 * @tparam F Type of functor
 * @tparam ArgsTupleT Type of tuple of arguments for functor
 *
 * @param functor Function to apply
 * @param err Error value to check against
 * @param args_tuple Tuple of arguments to pass to functor
 * @param log_result Whether the function result is already on the log-scale
 */
template <size_t ErrIndex, bool GradW7 = false, size_t NestedIndex = 0,
          typename F, typename ArgsTupleT>
double estimate_with_err_check(const F& functor, double err,
                               ArgsTupleT&& args_tuple,
                               bool log_result = true) {
  double result = math::apply([&](auto&&... args) { return functor(args...); },
                              args_tuple);
  double log_fabs_result = log_result ? log(fabs(result)) : fabs(result);
  if (log_fabs_result < err) {
    log_fabs_result = std::isinf(log_fabs_result) ? 0 : log_fabs_result;
    ArgsTupleT err_args_tuple = args_tuple;
    const double new_error
        = GradW7 ? err + log_fabs_result + LOG_TWO : err + log_fabs_result;
    assign_err<NestedIndex>(std::get<ErrIndex>(err_args_tuple), new_error);
    result = math::apply([&](auto&&... args) { return functor(args...); },
                         err_args_tuple);
  }
  return result;
}
}  // namespace internal

/**
 * Log-density function for the 5-parameter Wiener density.
 * See 'wiener_full_lpdf' for more comprehensive documentation
 *
 * @tparam T_y type of scalar
 * @tparam T_a type of boundary
 * @tparam T_t0 type of non-decision time
 * @tparam T_w type of relative starting point
 * @tparam T_v type of drift rate
 * @tparam T_sv type of inter-trial variability of drift rate
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param t0 The non-decision time
 * @param w The relative starting point
 * @param v The drift rate
 * @param sv The inter-trial variability of the drift rate
 * @param precision_derivatives Level of precision in estimation
 * @return The log of the Wiener first passage time density with
 *  the specified arguments for upper boundary responses
 */
template <bool propto = false, typename T_y, typename T_a, typename T_t0,
          typename T_w, typename T_v, typename T_sv>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv> wiener5_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const double& precision_derivatives) {
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_v_ref = ref_type_t<T_v>;
  using T_sv_ref = ref_type_t<T_sv>;

  const char* function_name = "wiener5_lpdf";
  if (size_zero(y, a, t0, w, v, sv)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v, T_sv>::value) {
    return 0;
  }

  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv);
  check_consistent_size(function_name, "Random variable", y, 1);
  check_consistent_size(function_name, "Boundary separation", a, 1);
  check_consistent_size(function_name, "Nondecision time", t0, 1);
  check_consistent_size(function_name, "A-priori bias", w, 1);
  check_consistent_size(function_name, "Drift rate", v, 1);
  check_consistent_size(function_name, "Inter-trial variability in drift rate",
                        sv, 1);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;
  T_sv_ref sv_ref = sv;

  check_positive_finite(function_name, "Random variable", value_of(y_ref));
  check_positive_finite(function_name, "Boundary separation", value_of(a_ref));
  check_nonnegative(function_name, "Nondecision time", value_of(t0_ref));
  check_finite(function_name, "Nondecision time", value_of(t0_ref));
  check_less(function_name, "A-priori bias", value_of(w_ref), 1);
  check_greater(function_name, "A-priori bias", value_of(w_ref), 0);
  check_finite(function_name, "Drift rate", value_of(v_ref));
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    value_of(sv_ref));
  check_finite(function_name, "Inter-trial variability in drift rate",
               value_of(sv_ref));

  const size_t N = max_size(y, a, t0, w, v, sv);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_a_ref> a_vec(a_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_w_ref> w_vec(w_ref);
  scalar_seq_view<T_v_ref> v_vec(v_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
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

  const double log_error_density = log(1e-6);
  const double log_error_derivative = log(precision_derivatives);
  const double log_error_absolute = log(1e-12);
  double density = 0.0;
  double log_density = 0.0;
  operands_and_partials<T_y_ref, T_a_ref, T_t0_ref, T_w_ref, T_v_ref, T_sv_ref>
      ops_partials(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    const double y_val = y_vec.val(i);
    const double a_val = a_vec.val(i);
    const double t0_val = t0_vec.val(i);
    const double w_val = w_vec.val(i);
    const double v_val = v_vec.val(i);
    const double sv_val = sv_vec.val(i);

    const auto params = std::make_tuple(y_val - t0_val, a_val, v_val, w_val,
                                        sv_val, log_error_absolute);

    density = internal::estimate_with_err_check<5>(
        [&](auto&&... args) { return internal::wiener5_density(args...); },
        log_error_density - LOG_TWO, params, false);
    log_density += density;

    const double new_est_err = density + log_error_derivative - LOG_FOUR;

    // computation of derivative for t and precision check in order to give
    // the value as deriv_y to edge1 and as -deriv_y to edge5
    const double deriv_y = internal::estimate_with_err_check<5>(
        [&](auto&&... args) { return internal::grad_wiener5_t(args...); },
        new_est_err, params);

    // computation of derivatives and precision checks
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[i] = deriv_y;
    }
    if (!is_constant_all<T_a>::value) {
      ops_partials.edge2_.partials_[i] = internal::estimate_with_err_check<5>(
          [&](auto&&... args) { return internal::grad_wiener5_a(args...); },
          new_est_err, params);
    }
    if (!is_constant_all<T_t0>::value) {
      ops_partials.edge3_.partials_[i] = -deriv_y;
    }
    if (!is_constant_all<T_w>::value) {
      ops_partials.edge4_.partials_[i] = internal::estimate_with_err_check<5>(
          [&](auto&&... args) { return internal::grad_wiener5_w(args...); },
          new_est_err, params);
    }
    if (!is_constant_all<T_v>::value) {
      ops_partials.edge5_.partials_[i] = internal::grad_wiener5_v(
          y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
    if (!is_constant_all<T_sv>::value) {
      ops_partials.edge6_.partials_[i] = internal::grad_wiener5_sv(
          y_val - t0_val, a_val, v_val, w_val, sv_val);
    }
  }  // end for loop
  return ops_partials.build(log_density);
}  // end wiener5_lpdf
}  // namespace math
}  // namespace stan
#endif
