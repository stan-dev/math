#ifndef STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER5_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename... Types>
using is_any_eigen_vector = disjunction<is_eigen_vector<Types>...>;
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
template <typename T_y, typename T_a, typename T_v,
          typename T_w, typename T_sv>
inline auto wiener5_compute_error_term(T_y&& y, T_a&& a, T_v&& v_value,
                                         T_w&& w_value, T_sv&& sv) noexcept {
  const auto w = 1.0 - w_value;
  const auto v = -v_value;
  const auto sv_sqr = square(sv);
  const auto one_plus_svsqr_y = 1 + sv_sqr * y;
  const auto two_avw = 2 * a * v * w;
  const auto two_log_a = 2 * log(a);
  return stan::math::eval((sv_sqr * square(a * w) - two_avw - square(v) * y)
                                / 2.0 / one_plus_svsqr_y
                            - two_log_a - 0.5 * log(one_plus_svsqr_y));
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
template <bool GradA, bool GradT, typename T_y, typename T_a,
          typename T_v_value, typename T_w_value, typename T_sv>
inline auto wiener5_density_part_one(T_y&& y, T_a&& a, T_v_value&& v_value,
                                       T_w_value&& w_value,
                                       T_sv&& sv) noexcept {
  const auto w = 1.0 - w_value;
  const auto v = -v_value;
  const auto sv_sqr = square(sv);
  const auto one_plus_svsqr_y = 1 + sv_sqr * y;
  if (GradT) {
      const auto two_avw = 2 * a * v * w;
      return -0.5
             * (square(sv_sqr) * (y + square(a * w)) + sv_sqr * (1 - two_avw)
                + square(v))
             / square(one_plus_svsqr_y);
  }
  const auto var_a = GradA ? w : a;
  const auto var_b = GradA ? a : w;
  return (-v * var_a + sv_sqr * square(var_a) * var_b) / one_plus_svsqr_y;
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
template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_small_t_plain(T_y&& y, T_a&& a, T_w_value&& w_value,
                                      T_err error) noexcept {
  const auto two_error = 2.0 * error;
  const auto y_asq = y / square(a);
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto w = 1.0 - w_value;

  auto u_eps =  fmin(-3.0,
                 (log(8.0) - log(27.0) + LOG_PI + 4.0 * log_y_asq + two_error));
  const auto arg = -3.0 * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  const auto n_1 = (sqrt(3.0 * y_asq) + w) / 2.0;
  const auto n_2 = (arg > 0) ? (0.5 * (sqrt(arg) - w)) : n_1;
  return ceil(fmax(n_1, n_2));
}

template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_small_t_gradw(T_y&& y, T_a&& a, T_w_value&& w_value,
                                      T_err error) noexcept {
  const auto two_error = 2.0 * error;
  const auto y_asq = y / square(a);
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto w = 1.0 - w_value;

  auto u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log_y_asq + two_error);
  const auto arg_mult = 1.0;
  const auto arg = -arg_mult * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));

  const auto n_1_factor = 3.0;
  const auto n_1 = (sqrt(n_1_factor * y_asq) + w) / 2.0;
  const auto n_2 = (arg > 0) ? (0.5 * (sqrt(arg) + w)) : n_1;
  return ceil(fmax(n_1, n_2));
}

template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_small_t_density(T_y&& y, T_a&& a, T_w_value&& w_value,
                                      T_err error) noexcept {
  const auto two_error = 2.0 * error;
  const auto y_asq = y / square(a);
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto w = 1.0 - w_value;

  auto u_eps = fmin(-1.0, LOG_TWO + LOG_PI + 2.0 * log_y_asq + two_error);
  const auto arg_mult = 1.0;
  const auto arg = -arg_mult * y_asq * (u_eps - sqrt(-2.0 * u_eps - 2.0));

  const auto n_1_factor = 2;
  const auto n_1 = (sqrt(n_1_factor * y_asq) + w) / 2.0;
  const auto n_2 = (arg > 0) ? (0.5 * (sqrt(arg) - w)) : n_1;
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
template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_largel_t_plain(T_y&& y, T_a&& a, T_w_value&& w_value,
                                       T_err error) noexcept {
  const auto y_asq = y / square(a);
  static constexpr double PI_SQUARED = pi() * pi();
  using scalar_t = return_type_t<T_y, T_a, T_w_value, T_err>;
  constexpr auto n_1_factor = 3.0;
  scalar_t n_1 = sqrt(n_1_factor / y_asq) / pi();
  const auto u_eps_arg
      = log(3.0) - log(5.0) + LOG_PI + 2.0 * log(y) - 2.0 * log(a) + error;
  const auto u_eps = fmin(-1.0, u_eps_arg);
  const auto arg_mult = (2.0 / PI_SQUARED / y_asq);
  const auto arg = -arg_mult * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  scalar_t n_2 = (arg > 0) ? sqrt(arg) : n_1;
  return ceil(fmax(n_1, n_2));

}

template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_largel_t_density(T_y&& y, T_a&& a, T_w_value&& w_value,
                                       T_err error) noexcept {
  const auto y_asq = y / square(a);
  static constexpr double PI_SQUARED = pi() * pi();
  using scalar_t = return_type_t<T_y, T_a, T_w_value, T_err>;
  scalar_t n_1 = 1.0 / (pi() * sqrt(y_asq));
  const auto two_log_piy = -2.0 * (LOG_PI + log(y) - 2.0 * log(a) + error);
  scalar_t n_2 = (two_log_piy >= 0) ? sqrt(two_log_piy / (PI_SQUARED * y_asq)) : 0.0;
  return ceil(fmax(n_1, n_2));
}

template <typename T_y, typename T_a,
          typename T_w_value, typename T_err>
inline auto wiener5_n_terms_largel_t_gradw(T_y&& y, T_a&& a, T_w_value&& w_value,
                                       T_err error) noexcept {
  const auto y_asq = y / square(a);
  static constexpr double PI_SQUARED = pi() * pi();
  using scalar_t = return_type_t<T_y, T_a, T_w_value, T_err>;
  constexpr auto n_1_factor = 2.0;
  scalar_t n_1 = sqrt(n_1_factor / y_asq) / pi();
  const auto u_eps_arg
      = log(4.0) - log(9.0) + 2.0 * LOG_PI + 3.0 * log(y) - 2.0 * log(a) + 2.0 * error;
  const auto u_eps = fmin(-1.0, u_eps_arg);
  const auto arg_mult = 1.0;
  const auto arg = -arg_mult * (u_eps - sqrt(-2.0 * u_eps - 2.0));
  scalar_t n_2 = (arg > 0) ? sqrt(arg / y_asq) / pi() : n_1;
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
template <typename T_y, typename T_a,
          typename T_w, typename T_nsmall, typename T_nlarge>
inline auto wiener5_log_sum_exp_gradw(T_y&& y, T_a&& a, T_w&& w_value,
                                T_nsmall&& n_terms_small_t,
                                T_nlarge&& n_terms_large_t) noexcept {
  const auto y_asq = y / square(a);
  reverse_pass_callback([w_value]{
      std::cout << std::endl;
      std::cout << "w_value partial log_sum val.val 2 = " << w_value.val().val() << std::endl;
      std::cout << "w_value partial log_sum val.adj 2 = " << w_value.val().adj() << std::endl;
      std::cout << "w_value partial log_sum d.val 2 = " << w_value.d().val() << std::endl;
      std::cout << "w_value partial log_sum d.adj 2 = " << w_value.d().adj() << std::endl;
      std::cout << std::endl;
  });
  const auto w = 1.0 - w_value;
  reverse_pass_callback([w, w_value]{
      std::cout << std::endl;
      std::cout << "w partial log_sum val.val 1 = " << w.val().val() << std::endl;
      std::cout << "w partial log_sum val.adj 1 = " << w.val().adj() << std::endl;
      std::cout << "w partial log_sum d.val 1 = " << w.d().val() << std::endl;
      std::cout << "w partial log_sum d.adj 1 = " << w.d().adj() << std::endl;
      std::cout << std::endl;
      std::cout << "w_value partial log_sum val.val 1 = " << w_value.val().val() << std::endl;
      std::cout << "w_value partial log_sum val.adj 1 = " << w_value.val().adj() << std::endl;
      std::cout << "w_value partial log_sum d.val 1 = " << w_value.d().val() << std::endl;
      std::cout << "w_value partial log_sum d.adj 1 = " << w_value.d().adj() << std::endl;
      std::cout << std::endl;
  });

  const bool small_n_terms_small_t = (2 * n_terms_small_t < n_terms_large_t);
  const auto scaling = small_n_terms_small_t ? inv(2.0 * y_asq) : y_asq / 2.0;
  using ret_t = return_type_t<T_y, T_a, T_w>;
  ret_t prev_val = NEGATIVE_INFTY;
  ret_t current_val = NEGATIVE_INFTY;
  int prev_sign = 1;
  int current_sign = 1;

  if (small_n_terms_small_t) {
    const auto mult = 3;
    const auto offset = y_asq;
    const auto sqrt_offset = sqrt(offset);
    for (auto k = n_terms_small_t; k >= 1; k--) {
      const auto wp2k = w + 2.0 * k;
      const auto wm2k = w - 2.0 * k;
      int wp2k_sign = (wp2k > sqrt_offset) ? 1 : -1;
      int wm2k_sign = (wm2k > sqrt_offset) ? 1 : -1;
      auto wp2k_quant = log(fabs((square(wp2k) - offset))) - square(wp2k) * scaling;
      auto wm2k_quant = log(fabs((square(wm2k) - offset))) - square(wm2k) * scaling;
      ret_t k_term;
      int k_sign;
      std::forward_as_tuple(k_term, k_sign) = log_sum_exp_signed(
          wm2k_quant, -1 * wm2k_sign, wp2k_quant, wp2k_sign);
      std::forward_as_tuple(current_val, current_sign)
          = log_sum_exp_signed(k_term, k_sign, prev_val, prev_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
    ret_t new_val = log(fabs(square(w) - offset)) - square(w) * scaling;
    reverse_pass_callback([w]{
      std::cout << "This Happened!!!" << std::endl;
      std::cout << "w partial log_sum val.val 3 = " << w.val().val() << std::endl;
      std::cout << "w partial log_sum val.adj 3 = " << w.val().adj() << std::endl;
      std::cout << "w partial log_sum d.val 3 = " << w.d().val() << std::endl;
      std::cout << "w partial log_sum d.adj 3 = " << w.d().adj() << std::endl;
      std::cout << std::endl;
    });
    int new_val_sign = (w > sqrt_offset ? 1 : -1);
    int factor_sign = 1;
    std::forward_as_tuple(current_val, current_sign)
        = log_sum_exp_signed(new_val, factor_sign * new_val_sign, current_val,
                             factor_sign * current_sign);
  } else {
    constexpr auto mult = 2.0;
    for (auto k = n_terms_large_t; k >= 1; k--) {
      const auto pi_k = k * pi();
      reverse_pass_callback([w]{
        std::cout << "w partial log_sum val.val 6 = " << w.val().val() << std::endl;
        std::cout << "w partial log_sum val.adj 6 = " << w.val().adj() << std::endl;
        std::cout << "w partial log_sum d.val 6 = " << w.d().val() << std::endl;
        std::cout << "w partial log_sum d.adj 6 = " << w.d().adj() << std::endl;
        std::cout << std::endl;
      });
      const auto check2 = pi_k * w;
      reverse_pass_callback([w, check2]{
        std::cout << "w partial log_sum val.val 5 = " << w.val().val() << std::endl;
        std::cout << "w partial log_sum val.adj 5 = " << w.val().adj() << std::endl;
        std::cout << "w partial log_sum d.val 5 = " << w.d().val() << std::endl;
        std::cout << "w partial log_sum d.adj 5 = " << w.d().adj() << std::endl;
        std::cout << std::endl;
        std::cout << "check2 partial log_sum val.val 5 = " << check2.val().val() << std::endl;
        std::cout << "check2 partial log_sum val.adj 5 = " << check2.val().adj() << std::endl;
        std::cout << "check2 partial log_sum d.val 5 = " << check2.d().val() << std::endl;
        std::cout << "check2 partial log_sum d.adj 5 = " << check2.d().adj() << std::endl;
        std::cout << std::endl;
      });
      const auto check = cos(check2);
      reverse_pass_callback([w, check]{
        std::cout << "w partial log_sum val.val 4 = " << w.val().val() << std::endl;
        std::cout << "w partial log_sum val.adj 4 = " << w.val().adj() << std::endl;
        std::cout << "w partial log_sum d.val 4 = " << w.d().val() << std::endl;
        std::cout << "w partial log_sum d.adj 4 = " << w.d().adj() << std::endl;
        std::cout << std::endl;
        std::cout << "check partial log_sum val.val 5 = " << check.val().val() << std::endl;
        std::cout << "check partial log_sum val.adj 5 = " << check.val().adj() << std::endl;
        std::cout << "check partial log_sum d.val 5 = " << check.d().val() << std::endl;
        std::cout << "check partial log_sum d.adj 5 = " << check.d().adj() << std::endl;
        std::cout << std::endl;
      });
      int check_sign = sign(check);
      reverse_pass_callback([check]{
        std::cout << "check partial log_sum val.val 7 = " << check.val().val() << std::endl;
        std::cout << "check partial log_sum val.adj 7 = " << check.val().adj() << std::endl;
        std::cout << "check partial log_sum d.val 7 = " << check.d().val() << std::endl;
        std::cout << "check partial log_sum d.adj 7 = " << check.d().adj() << std::endl;
        std::cout << std::endl;
      });
      auto check3 = fabs(check);
      reverse_pass_callback([check3]{
        std::cout << "check3 partial log_sum val.val 8 = " << check3.val().val() << std::endl;
        std::cout << "check3 partial log_sum val.adj 8 = " << check3.val().adj() << std::endl;
        std::cout << "check3 partial log_sum d.val 8 = " << check3.d().val() << std::endl;
        std::cout << "check3 partial log_sum d.adj 8 = " << check3.d().adj() << std::endl;
        std::cout << std::endl;
      });
      auto check4 = log(check3);
      reverse_pass_callback([k, check, check4] {
        std::cout << "k: " << k << std::endl;
        std::cout << "check4 partial log_sum val.val 6 = " << check4.val().val() << std::endl;
        std::cout << "check4 partial log_sum val.adj 6 = " << check4.val().adj() << std::endl;
        std::cout << "check4 partial log_sum d.val 6 = " << check4.d().val() << std::endl;
        std::cout << "check4 partial log_sum d.adj 6 = " << check4.d().adj() << std::endl;
        std::cout << std::endl;
      });
      auto n_terms_large_t_quant
          = mult * log(k) - square(pi_k) * scaling + check4;
      reverse_pass_callback([k, check, check3] {
        std::cout << "k: " << k << std::endl;
        std::cout << "check3 partial log_sum val.val 9 = " << check3.val().val() << std::endl;
        std::cout << "check3 partial log_sum val.adj 9 = " << check3.val().adj() << std::endl;
        std::cout << "check3 partial log_sum d.val 9 = " << check3.d().val() << std::endl;
        std::cout << "check3 partial log_sum d.adj 9 = " << check3.d().adj() << std::endl;
        std::cout << std::endl;
        std::cout << "check partial log_sum val.val 6 = " << check.val().val() << std::endl;
        std::cout << "check partial log_sum val.adj 6 = " << check.val().adj() << std::endl;
        std::cout << "check partial log_sum d.val 6 = " << check.d().val() << std::endl;
        std::cout << "check partial log_sum d.adj 6 = " << check.d().adj() << std::endl;
        std::cout << std::endl;
      });
      std::forward_as_tuple(current_val, current_sign) = log_sum_exp_signed(
          prev_val, prev_sign, n_terms_large_t_quant, check_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
  }
  return std::make_pair(current_val, current_sign);
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
template <typename T_y, typename T_a,
          typename T_w, typename T_nsmall, typename T_nlarge>
inline auto wiener5_log_sum_exp_density(T_y&& y, T_a&& a, T_w&& w_value,
                                T_nsmall&& n_terms_small_t,
                                T_nlarge&& n_terms_large_t) noexcept {
  const auto y_asq = y / square(a);
  const auto w = 1.0 - w_value;
  const bool small_n_terms_small_t = (2 * n_terms_small_t <= n_terms_large_t);
  const auto scaling = small_n_terms_small_t ? inv(2.0 * y_asq) : y_asq / 2.0;
  using ret_t = return_type_t<T_y, T_a, T_w>;
  ret_t prev_val = NEGATIVE_INFTY;
  ret_t current_val = NEGATIVE_INFTY;
  int prev_sign = 1;
  int current_sign = 1;

  if (small_n_terms_small_t) {
    const auto mult = 1;
    const auto offset = 0;
    const auto sqrt_offset = sqrt(offset);
    for (auto k = n_terms_small_t; k >= 1; k--) {
      const auto wp2k = w + 2.0 * k;
      const auto wm2k = w - 2.0 * k;
      int wp2k_sign = (wp2k > sqrt_offset) ? 1 : -1;
      int wm2k_sign = (wm2k > sqrt_offset) ? 1 : -1;
      auto wp2k_quant = mult * log(wp2k_sign * wp2k) - square(wp2k) * scaling;
      auto wm2k_quant = mult * log(wm2k_sign * wm2k) - square(wm2k) * scaling;
      ret_t k_term;
      int k_sign;
      std::forward_as_tuple(k_term, k_sign) = log_sum_exp_signed(
          wm2k_quant, -1 * wm2k_sign, wp2k_quant, wp2k_sign);
      std::forward_as_tuple(current_val, current_sign)
          = log_sum_exp_signed(k_term, k_sign, prev_val, prev_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
    ret_t new_val = mult * log(w) - square(w) * scaling;
    int new_val_sign = (new_val > 0 ? 1 : -1);
    constexpr int factor_sign = -1;
    std::forward_as_tuple(current_val, current_sign)
        = log_sum_exp_signed(new_val, factor_sign * new_val_sign, current_val,
                             factor_sign * current_sign);
  } else {
    auto mult = 1.0;
    for (auto k = n_terms_large_t; k >= 1; k--) {
      const auto pi_k = k * pi();
      const auto check = sin(pi_k * w);
      int check_sign = sign(check);
      auto n_terms_large_t_quant
          = mult * log(k) - square(pi_k) * scaling + log(fabs(check));
      std::forward_as_tuple(current_val, current_sign) = log_sum_exp_signed(
          prev_val, prev_sign, n_terms_large_t_quant, check_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
  }
  return std::make_pair(current_val, current_sign);
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
template <typename T_y, typename T_a,
          typename T_w, typename T_nsmall, typename T_nlarge>
inline auto wiener5_log_sum_exp_plain(T_y&& y, T_a&& a, T_w&& w_value,
                                T_nsmall&& n_terms_small_t,
                                T_nlarge&& n_terms_large_t) noexcept {
  const auto y_asq = y / square(a);
  const auto w = 1.0 - w_value;
  const bool small_n_terms_small_t = (2 * n_terms_small_t < n_terms_large_t);
  const auto scaling = small_n_terms_small_t ? inv(2.0 * y_asq) : y_asq / 2.0;
  using ret_t = return_type_t<T_y, T_a, T_w>;
  ret_t prev_val = NEGATIVE_INFTY;
  ret_t current_val = NEGATIVE_INFTY;
  int prev_sign = 1;
  int current_sign = 1;

  if (small_n_terms_small_t) {
    constexpr auto mult = 3.0;
    constexpr auto offset = 0.0;
    const auto sqrt_offset = sqrt(offset);
    for (auto k = n_terms_small_t; k >= 1; k--) {
      const auto wp2k = w + 2.0 * k;
      const auto wm2k = w - 2.0 * k;
      int wp2k_sign = (wp2k > sqrt_offset) ? 1 : -1;
      int wm2k_sign = (wm2k > sqrt_offset) ? 1 : -1;
      auto wp2k_quant = mult * log(wp2k_sign * wp2k) - square(wp2k) * scaling;
      auto wm2k_quant = log(wm2k_sign * wm2k) - square(wm2k) * scaling;
      ret_t k_term;
      int k_sign;
      std::forward_as_tuple(k_term, k_sign) = log_sum_exp_signed(
          wm2k_quant, -1 * wm2k_sign, wp2k_quant, wp2k_sign);
      std::forward_as_tuple(current_val, current_sign)
          = log_sum_exp_signed(k_term, k_sign, prev_val, prev_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
    ret_t new_val = mult * log(w) - square(w) * scaling;
    int new_val_sign = (new_val > 0 ? 1 : -1);
    constexpr int factor_sign = -1;
    std::forward_as_tuple(current_val, current_sign)
        = log_sum_exp_signed(new_val, factor_sign * new_val_sign, current_val,
                             factor_sign * current_sign);
  } else {
    auto mult = 3.0;
    for (auto k = n_terms_large_t; k >= 1; k--) {
      const auto pi_k = k * pi();
      const auto check = sin(pi_k * w);
      int check_sign = sign(check);
      auto n_terms_large_t_quant
          = mult * log(k) - square(pi_k) * scaling + log(fabs(check));
      std::forward_as_tuple(current_val, current_sign) = log_sum_exp_signed(
          prev_val, prev_sign, n_terms_large_t_quant, check_sign);
      prev_val = current_val;
      prev_sign = current_sign;
    }
  }
  return std::make_pair(current_val, current_sign);
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
template <bool NaturalScale = false, typename T_y,
          typename T_a, typename T_w, typename T_v, typename T_sv, typename T_err>
inline auto wiener5_density(const T_y& y, const T_a& a, const T_v& v_value,
                              const T_w& w_value, const T_sv& sv,
                              T_err err = log(1e-12)) noexcept {
  const auto error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const auto error = (err - error_term);
  const auto n_terms_small_t
      = wiener5_n_terms_small_t_density(y, a, w_value, error);
  const auto n_terms_large_t
      = wiener5_n_terms_largel_t_density(y, a, w_value, error);

  auto res = wiener5_log_sum_exp_density(
                 y, a, w_value, n_terms_small_t, n_terms_large_t)
                 .first;
  if (2 * n_terms_small_t <= n_terms_large_t) {
    auto log_density = error_term - 0.5 * LOG_TWO - LOG_SQRT_PI
                         - 1.5 * (log(y) - 2 * log(a)) + res;
    return NaturalScale ? exp(log_density) : log_density;
  } else {
    auto log_density = error_term + res + LOG_PI;
    return NaturalScale ? exp(log_density) : log_density;
  }
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 't'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 * @tparam Scalar type of scalars
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. t
 */
template <bool WrtLog = false, typename T_y, typename T_a,
          typename T_w, typename T_v, typename T_sv, typename T_err>
inline auto wiener5_grad_t(const T_y& y, const T_a& a, const T_v& v_value,
                             const T_w& w_value, const T_sv& sv,
                             T_err err = log(1e-12)) noexcept {
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const auto density_part_one = wiener5_density_part_one<false, true>(
      y, a, v_value, w_value, sv);
  const auto error = (err - error_term) + two_log_a;

  const auto n_terms_small_t
      = wiener5_n_terms_small_t_plain(y, a, w_value, error);
  const auto n_terms_large_t
      = wiener5_n_terms_largel_t_plain(y, a, w_value, error);
  auto grad_result = wiener5_log_sum_exp_plain(
          y, a, w_value, n_terms_small_t, n_terms_large_t);
  auto result = grad_result.first;
  int newsign = grad_result.second;

  const auto error_log_density
      = log(fmax(fabs(density_part_one - 1.5 / y), fabs(density_part_one)));
  const auto log_density = wiener5_density<false>(
      y, a, v_value, w_value, sv, err - error_log_density);
  if (2 * n_terms_small_t < n_terms_large_t) {
    const auto ans = density_part_one - 1.5 / y
          + newsign
                * exp(error_term - two_log_a - 1.5 * LOG_TWO - LOG_SQRT_PI
                      - 3.5 * log_y_asq + result - log_density);
  return WrtLog ? ans * exp(log_density) : ans;
  } else {
    const auto ans = density_part_one
          - newsign
                * exp(error_term - two_log_a + 3.0 * LOG_PI - LOG_TWO + result
                      - log_density);
  return WrtLog ? ans * exp(log_density) : ans;
  }
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'a'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 * @tparam Scalar type of scalars
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. a
 */
template <bool WrtLog = false, typename Scalar, typename T_y, typename T_a,
          typename T_w, typename T_v, typename T_sv>
inline Scalar wiener5_grad_a(const T_y& y, const T_a& a, const T_v& v_value,
                             const T_w& w_value, const T_sv& sv,
                             Scalar err = log(1e-12)) noexcept {
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const auto density_part_one = wiener5_density_part_one<true, false>(
      y, a, v_value, w_value, sv);
  const auto error = err - error_term + 3 * log(a) - log(y) - LOG_TWO;

  const auto n_terms_small_t
      = wiener5_n_terms_small_t_plain(y, a, w_value, error);
  const auto n_terms_large_t
      = wiener5_n_terms_largel_t_plain(y, a, w_value, error);
  Scalar result;
  int newsign;
  std::forward_as_tuple(result, newsign)
      = wiener5_log_sum_exp_plain(
          y, a, w_value, n_terms_small_t, n_terms_large_t);

  const auto error_log_density = log(
      fmax(fabs(density_part_one + 1.0 / a), fabs(density_part_one - 2.0 / a)));
  const auto log_density = wiener5_density<false>(
      y, a, v_value, w_value, sv, err - error_log_density);
  Scalar ans;
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
 * @tparam Scalar type of scalars
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. v
 */
template <bool WrtLog = false, typename Scalar, typename T_y, typename T_a,
          typename T_w, typename T_v, typename T_sv>
inline Scalar wiener5_grad_v(const T_y& y, const T_a& a, const T_v& v_value,
                             const T_w& w_value, const T_sv& sv,
                             Scalar err = log(1e-12)) noexcept {
  Scalar ans = (a * (1 - w_value) - v_value * y) / (1 + square(sv) * y);
  return WrtLog ? ans
                      * wiener5_density<true>(y, a, v_value, w_value,
                                                      sv, err)
                : ans;
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'w'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 * @tparam Scalar type of scalars
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. w
 */
template <bool WrtLog = false, typename T_y, typename T_a,
          typename T_w, typename T_v, typename T_sv, typename T_err>
inline auto wiener5_grad_w(const T_y& y, const T_a& a, const T_v& v_value,
                             const T_w& w_value, const T_sv& sv,
                             T_err err = log(1e-12)) noexcept {
  const auto two_log_a = 2 * log(a);
  const auto log_y_asq = log(y) - two_log_a;
  const auto error_term
      = wiener5_compute_error_term(y, a, v_value, w_value, sv);
  const auto density_part_one
      = wiener5_density_part_one<false, false>(y, a, v_value, w_value,
                                                       sv);
  reverse_pass_callback([w_value]{
      std::cout << "w partial 6 = " << w_value.val().adj() << std::endl;
  });
  const auto error = (err - error_term);

  reverse_pass_callback([w_value]{
      std::cout << "w partial 5 = " << w_value.val().adj() << std::endl;
  });
  const auto n_terms_small_t
      = wiener5_n_terms_small_t_gradw(y, a, w_value, error);
  reverse_pass_callback([w_value]{
      std::cout << "w partial 4 = " << w_value.val().adj() << std::endl;
  });
  const auto n_terms_large_t
      = wiener5_n_terms_largel_t_gradw(y, a, w_value, error);
  reverse_pass_callback([w_value]{
      std::cout << "w partial 3 = " << w_value.val().adj() << std::endl;
  });
  auto grad_result
      = wiener5_log_sum_exp_gradw(y, a, w_value, n_terms_small_t,
                                                 n_terms_large_t);
  auto result = grad_result.first;
  auto newsign = grad_result.second;
  reverse_pass_callback([w_value]{
      std::cout << "w partial 2 = " << w_value.val().adj() << std::endl;
  });
  const auto log_density = wiener5_density<false>(
      y, a, v_value, w_value, sv, err - log(fabs(density_part_one)));
  reverse_pass_callback([w_value]{
      std::cout << "w partial 1 = " << w_value.val().adj() << std::endl;
  });
  if (2 * n_terms_small_t < n_terms_large_t) {
    const auto ans = -(density_part_one
            - newsign
                  * exp(result - (log_density - error_term) - 2.5 * log_y_asq
                        - 0.5 * LOG_TWO - 0.5 * LOG_PI));
  return WrtLog ? ans * exp(log_density) : ans;
  } else {
    const auto ans = -(density_part_one
            + newsign * exp(result - (log_density - error_term) + 2 * LOG_PI));
  return WrtLog ? ans * exp(log_density) : ans;
  }
}

/**
 * Calculate the derivative of the wiener5 density w.r.t. 'sv'
 *
 * @tparam WrtLog Whether to return the derivative w.r.t.
 *                  the natural (true) or log-scale (false) density
 * @tparam Scalar type of scalars
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v_value The drift rate
 * @param w_value The relative starting point
 * @param sv The inter-trial variability of the drift rate
 * @param err The log error tolerance
 * @return Gradient w.r.t. sv
 */
template <bool WrtLog = false, typename Scalar, typename T_y, typename T_a,
          typename T_w, typename T_v, typename T_sv>
inline Scalar wiener5_grad_sv(const T_y& y, const T_a& a, const T_v& v_value,
                              const T_w& w_value, const T_sv& sv,
                              Scalar err = log(1e-12)) noexcept {
  const auto one_plus_svsqr_y = 1 + square(sv) * y;
  const auto w = 1.0 - w_value;
  const auto v = -v_value;
  const auto t1 = -y / one_plus_svsqr_y;
  const auto t2 = (square(a * w) + 2 * a * v * w * y + square(v * y))
                    / square(one_plus_svsqr_y);
  const auto ans = sv * (t1 + t2);
  return WrtLog ? ans
                      * wiener5_density<true>(y, a, v_value, w_value,
                                                      sv, err)
                : ans;
}

/**
 * Utility function for replacing a value with a specified error value
 */
template <size_t NestedIndex, typename Scalar1, typename Scalar2>
inline void assign_err(Scalar1 arg, Scalar2 err) {
  arg = err;
}

/**
 * Utility function for replacing a value with a specified error value,
 * overload for use when the value is stored within a tuple.
 */
template <size_t NestedIndex, typename Scalar, typename... TArgs>
inline void assign_err(std::tuple<TArgs...>& args_tuple, Scalar err) {
  std::get<NestedIndex>(args_tuple) = err;
}

/**
 * Utility function for estimating a function with a given set of arguments,
 * checking the result against a provided error tolerance, and re-estimating
 * the function with the increased error if it fails.
 *
 * @tparam Scalar type of scalars
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
template <size_t ErrIndex, bool GradW7 = false,
          size_t NestedIndex = 0, bool LogResult = true, typename F, typename T_err,
          typename... ArgsTupleT>
inline auto estimate_with_err_check(const F& functor, T_err err,
                               ArgsTupleT&&... args_tuple) {
  auto result = functor(args_tuple...);
  auto log_fabs_result = LogResult ? log(fabs(result)) : fabs(result);
  if (log_fabs_result < err) {
    log_fabs_result = is_inf(log_fabs_result) ? 0 : log_fabs_result;
    auto err_args_tuple = std::make_tuple(args_tuple...);
    const auto new_error
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
          typename T_w, typename T_v, typename T_sv,
          typename ReturnT = return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv>>
inline return_type_t<T_y, T_a, T_t0, T_w, T_v, T_sv> wiener5_lpdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
    const T_sv& sv, const double& precision_derivatives) {
  using T_partials_return = partials_type_t<ReturnT>;

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v, T_sv>::value) {
    return 0;
  }

  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_a_ref = ref_type_if_t<!is_constant<T_a>::value, T_a>;
  using T_t0_ref = ref_type_if_t<!is_constant<T_t0>::value, T_t0>;
  using T_w_ref = ref_type_if_t<!is_constant<T_w>::value, T_w>;
  using T_v_ref = ref_type_if_t<!is_constant<T_v>::value, T_v>;
  using T_sv_ref = ref_type_if_t<!is_constant<T_sv>::value, T_sv>;

  static constexpr const char* function_name = "wiener5_lpdf";

  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv);

  T_y_ref y_ref = y;
  T_a_ref a_ref = a;
  T_t0_ref t0_ref = t0;
  T_w_ref w_ref = w;
  T_v_ref v_ref = v;
  T_sv_ref sv_ref = sv;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) a_val = to_ref(as_value_column_array_or_scalar(a_ref));
  decltype(auto) v_val = to_ref(as_value_column_array_or_scalar(v_ref));
  decltype(auto) w_val = to_ref(as_value_column_array_or_scalar(w_ref));
  decltype(auto) t0_val = to_ref(as_value_column_array_or_scalar(t0_ref));
  decltype(auto) sv_val = to_ref(as_value_column_array_or_scalar(sv_ref));
  check_positive_finite(function_name, "Random variable", y_val);
  check_positive_finite(function_name, "Boundary separation", a_val);
  check_finite(function_name, "Drift rate", v_val);
  check_less(function_name, "A-priori bias", w_val, 1);
  check_greater(function_name, "A-priori bias", w_val, 0);
  check_nonnegative(function_name, "Nondecision time", t0_val);
  check_finite(function_name, "Nondecision time", t0_val);
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    sv_val);
  check_finite(function_name, "Inter-trial variability in drift rate", sv_val);

  if (size_zero(y, a, t0, w, v, sv)) {
    return 0;
  }
  const size_t N = max_size(y, a, t0, w, v, sv);
  if (!N) {
    return 0;
  }

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
  const partials_type_t<ReturnT> log_error_absolute = log(1e-12);
  T_partials_return log_density = 0.0;
  auto ops_partials
      = make_partials_propagator(y_ref, a_ref, t0_ref, w_ref, v_ref, sv_ref);

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    const auto y_value = y_vec.val(i);
    const auto a_value = a_vec.val(i);
    const auto t0_value = t0_vec.val(i);
    const auto w_value = w_vec.val(i);
    const auto v_value = v_vec.val(i);
    const auto sv_value = sv_vec.val(i);
    reverse_pass_callback([w_value]{
      std::cout << "w end loop = " << w_value.val().adj() << std::endl;
    });

    auto l_density = internal::estimate_with_err_check<5, false, 0, false>(
        [&](auto&&... args) {
          return internal::wiener5_density<false>(args...);
        },
        log_error_density - LOG_TWO, y_value - t0_value, a_value, v_value,
        w_value, sv_value, log_error_absolute);

    log_density += l_density;

    const auto new_est_err = l_density + log_error_derivative - LOG_FOUR;

    if (!is_constant_all<T_y, T_t0>::value) {
      // computation of derivative for t and precision check in order to give
      // the value as deriv_y to edge1 and as -deriv_y to edge5
      const auto deriv_y = internal::estimate_with_err_check<5, false, 0, true>(
          [&](auto&&... args) {
            return internal::wiener5_grad_t<false>(args...);
          },
          new_est_err, y_value - t0_value, a_value, v_value, w_value, sv_value,
          log_error_absolute);
      // computation of derivatives and precision checks
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[i] = deriv_y;
      }
      if (!is_constant_all<T_t0>::value) {
        partials<2>(ops_partials)[i] = -deriv_y;
      }
    }
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i] = internal::estimate_with_err_check<5, false, 0, true>(
          [&](auto&&... args) {
            return internal::wiener5_grad_a<false, T_partials_return>(args...);
          },
          new_est_err, y_value - t0_value, a_value, v_value, w_value, sv_value,
          log_error_absolute);
    }
    reverse_pass_callback([w_value]{
      std::cout << "w partial = " << w_value.val().adj() << std::endl;
    });
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i] = internal::estimate_with_err_check<5, false, 0, true>(
          [&](auto&&... args) {
            return internal::wiener5_grad_w<false>(args...);
          },
          new_est_err, y_value - t0_value, a_value, v_value, w_value, sv_value,
          log_error_absolute);
    }
    reverse_pass_callback([w_value]{
      std::cout << "w pre partial = " << w_value.val().adj() << std::endl;
    });
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener5_grad_v<false, T_partials_return>(
              y_value - t0_value, a_value, v_value, w_value, sv_value);
    }
    if (!is_constant_all<T_sv>::value) {
      partials<5>(ops_partials)[i]
          = internal::wiener5_grad_sv<false, T_partials_return>(
              y_value - t0_value, a_value, v_value, w_value, sv_value);
    }
  }  // end for loop
  return ops_partials.build(log_density);
}  // end wiener5_lpdf
}  // namespace math
}  // namespace stan
#endif
