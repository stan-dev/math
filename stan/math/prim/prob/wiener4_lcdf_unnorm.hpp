#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCDF_UNNORM_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCDF_UNNORM_HPP

#include <stan/math/prim/prob/wiener5_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Make the expression finite
 *
 * @param x Expression to test
 * @return Expression or limited to maximum numeric_limit
 */
template <typename T_x>
inline auto make_finite(const T_x& x) {
  if (x < std::numeric_limits<T_x>::max()) {
    return x;
  } else {
    return std::numeric_limits<T_x>::max();
  }
}

/**
 * Calculate the probability term 'P' on log scale for distribution
 *
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @return 'P' term
 */
template <typename T_a, typename T_w, typename T_v>
inline auto log_probability_distribution(const T_a& a, const T_v& v,
                                         const T_w& w) {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  if (fabs(v) == 0.0) {
    return ret_t(log1m(w));
  }
  auto two_va = 2.0 * v * a;
  auto minus_two_va_one_minus_w = -two_va * (1.0 - w);
  // This split prevents abort errors
  if (minus_two_va_one_minus_w < 0) {
    const auto exp_arg = exp(minus_two_va_one_minus_w);
    auto two_vaw = two_va * w;
    if (two_vaw > minus_two_va_one_minus_w) {
      return log1m(exp_arg) - log_diff_exp(two_vaw, minus_two_va_one_minus_w);
    } else if (two_vaw < minus_two_va_one_minus_w) {
      return log1m(exp_arg) - log_diff_exp(minus_two_va_one_minus_w, two_vaw);
    } else {
      return log1m(exp_arg) - NEGATIVE_INFTY;
    }
  } else {
    return log1m_exp(-minus_two_va_one_minus_w) - log1m_exp(two_va);
  }
}

/**
 * Calculate the probability term 'P' on log scale for grad_a and grad_v
 *
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @return 'P' term
 */
template <typename T_a, typename T_w, typename T_v>
inline auto log_probability_GradAV(const T_a& a, const T_v& v, const T_w& w) {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  if (fabs(v) == 0.0) {
    return ret_t(-w);
  }
  auto nearly_one = ret_t(1.0 - 1.1 * 1.0e-5);
  ret_t log_prob;
  if (v < 0) {
    const auto two_av = 2.0 * a * v;
    const auto two_va_one_minus_w = (two_av * (1.0 - w));
    const auto two_avw = two_av * w;
    const auto exp_two_va_one_minus_w = exp(two_va_one_minus_w);
    const auto exp_two_avw = exp(two_avw);
    const auto exp_two_av = exp(two_av);
    if (((exp_two_va_one_minus_w >= nearly_one) || (exp_two_avw >= nearly_one))
        || (exp_two_av >= nearly_one)) {
      return ret_t(-w);
    }
    log_prob = LOG_TWO + two_va_one_minus_w - log1m(exp_two_va_one_minus_w);
    auto log_quotient = log1m(exp_two_avw) - log1m(exp_two_av);
    if (log(w) > log_quotient) {
      return exp(log_prob) * (w - exp(log_quotient));
    } else {
      return -exp(log_prob) * (exp(log_quotient) - w);
    }
  } else {
    const auto minus_two_av = -2.0 * a * v;
    const auto minus_two_va_one_minus_w = minus_two_av * (1.0 - w);
    const auto exp_minus_two_va_one_minus_w = exp(minus_two_va_one_minus_w);
    const auto exp_minus_two_av = exp(minus_two_av);
    if ((exp_minus_two_va_one_minus_w >= nearly_one)
        || (exp_minus_two_av >= nearly_one)) {
      return ret_t(-w);
    }
    log_prob = LOG_TWO - log1m(exp_minus_two_va_one_minus_w);
    ret_t log_quotient;
    if (minus_two_va_one_minus_w > minus_two_av) {
      log_quotient = log_diff_exp(minus_two_va_one_minus_w, minus_two_av)
                     - log1m(exp_minus_two_av);
    } else if (minus_two_va_one_minus_w < minus_two_av) {
      log_quotient = log_diff_exp(minus_two_av, minus_two_va_one_minus_w)
                     - log1m(exp_minus_two_av);
    } else {
      log_quotient = NEGATIVE_INFTY;
    }
    if (log(w) > log_quotient) {
      return -exp(log_prob + log_diff_exp(log(w), log_quotient));
    } else {
      return exp(log_prob + log_diff_exp(log_quotient, log(w)));
    }
  }
}

/**
 * Log of Mill's ratio for the normal distribution
 *
 * @param x A scalar
 * @return The natural logarithm of Mill's ratio
 */
template <typename T_x>
inline auto logMill(T_x&& x) {
  return std_normal_lcdf(-x) - std_normal_lpdf(x);
}

/**
 * Calculate the wiener4 distribution
 *
 * @tparam NaturalScale Whether to return the distribution on natural or
 * log-scale
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return distribution
 */
template <bool NaturalScale = false, typename T_y, typename T_a, typename T_w,
          typename T_v, typename T_err>
inline auto wiener4_distribution(const T_y& y, const T_a& a, const T_v& v,
                                 const T_w& w, T_err log_err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto neg_v = -v;
  const auto one_m_w = 1.0 - w;

  const auto one_m_w_a_neg_v = one_m_w * a * neg_v;

  const auto K1 = 0.5 * (fabs(neg_v) / a * y - one_m_w);
  const auto arg = fmax(
      0.0, fmin(1.0, exp(one_m_w_a_neg_v + square(neg_v) * y / 2.0 + log_err)
                         / 2.0));
  const auto K2 = (arg == 0) ? INFTY
                             : (arg == 1) ? NEGATIVE_INFTY
                                          : -sqrt(y) / 2.0 / a * inv_Phi(arg);
  const auto K_small_value = ceil(fmax(K1, K1 + K2));

  const auto api = a / pi();
  const auto v_square = square(neg_v);
  const auto sqrtL1 = sqrt(1.0 / y) * api;
  const auto sqrtL2 = sqrt(fmax(
      1.0, -2.0 / y * square(api)
               * (log_err + log(pi() * y / 2.0 * (v_square + square(pi() / a)))
                  + one_m_w_a_neg_v + v_square * y / 2.0)));
  const auto K_large_value = ceil(fmax(sqrtL1, sqrtL2));

  auto lg = LOG_TWO + LOG_PI - 2.0 * log(a);

  if (3 * K_small_value < K_large_value) {
    const auto sqrt_y = sqrt(y);
    const auto neg_vy = neg_v * y;
    ret_t fplus = NEGATIVE_INFTY;
    ret_t fminus = NEGATIVE_INFTY;
    for (auto k = K_small_value; k >= 0; --k) {
      auto rj = a * (2.0 * k + one_m_w);
      auto dj = std_normal_lpdf(rj / sqrt_y);
      auto pos1 = dj + logMill((rj - neg_vy) / sqrt_y);
      auto pos2 = dj + logMill((rj + neg_vy) / sqrt_y);
      fplus = log_sum_exp(fplus, log_sum_exp(pos1, pos2));
      rj = a * (2.0 * k + 2.0 - one_m_w);
      dj = std_normal_lpdf(rj / sqrt_y);
      auto neg1 = dj + logMill((rj - neg_vy) / sqrt_y);
      auto neg2 = dj + logMill((rj + neg_vy) / sqrt_y);
      fminus = log_sum_exp(fminus, log_sum_exp(neg1, neg2));
    }
    auto ans = ret_t(0.0);
    ans = fplus > fminus ? log_diff_exp(fplus, fminus)
                         : log_diff_exp(fminus, fplus);
    ret_t log_distribution = ans - one_m_w_a_neg_v - square(neg_v) * y / 2;
    return NaturalScale ? exp(log_distribution) : log_distribution;
  }
  const auto log_a = log(a);
  const auto log_v = log(fabs(neg_v));
  ret_t fplus = NEGATIVE_INFTY;
  ret_t fminus = NEGATIVE_INFTY;
  for (auto k = K_large_value; k > 0; --k) {
    auto log_k = log(k);
    auto k_pi = k * pi();
    auto sin_k_pi_w = sin(k_pi * one_m_w);
    if (sin_k_pi_w > 0) {
      fplus = log_sum_exp(
          fplus, log_k
                     - log_sum_exp(2.0 * log_v, 2.0 * (log_k + LOG_PI - log_a))
                     - 0.5 * square(k_pi / a) * y + log(sin_k_pi_w));
    } else if (sin_k_pi_w < 0) {
      fminus = log_sum_exp(
          fminus, log_k
                      - log_sum_exp(2.0 * log_v, 2.0 * (log_k + LOG_PI - log_a))
                      - 0.5 * square(k_pi / a) * y + log(-sin_k_pi_w));
    }
  }
  ret_t ans = NEGATIVE_INFTY;
  ans = fplus > fminus ? log_diff_exp(fplus, fminus)
                       : log_diff_exp(fminus, fplus);
  auto summand_1 = log_probability_distribution(a, neg_v, one_m_w);
  auto summand_2 = lg + (ans - one_m_w_a_neg_v - 0.5 * square(neg_v) * y);
  ret_t log_distribution = NEGATIVE_INFTY;
  if (summand_1 > summand_2) {
    log_distribution = log_diff_exp(summand_1, summand_2);
  } else if (summand_1 < summand_2) {
    log_distribution = log_diff_exp(summand_2, summand_1);
  }
  return NaturalScale ? exp(log_distribution) : log_distribution;
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'a' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The value of the distribution
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to a
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_a(const T_y& y, const T_a& a, const T_v& v,
                               const T_w& w, T_cdf&& cdf,
                               T_err log_err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto neg_v = -v;
  const auto one_m_w = 1 - w;

  const auto one_m_w_neg_v = one_m_w * neg_v;
  const auto one_m_w_a_neg_v = one_m_w_neg_v * a;

  const auto log_y = log(y);
  const auto log_a = log(a);
  auto C1 = ret_t(
      LOG_TWO - log_sum_exp(2.0 * log(fabs(neg_v)), 2.0 * (LOG_PI - log_a)));
  C1 = log_sum_exp(C1, log_y);
  const auto factor = one_m_w_a_neg_v + square(neg_v) * y / 2.0 + log_err;
  const auto alphK = fmin(factor + LOG_PI + log_y + log_a - LOG_TWO - C1, 0.0);
  const auto K = a / pi() / sqrt(y);
  const auto K_large_value
      = ceil(fmax(fmax(sqrt(-2.0 * alphK / y) * a / pi(), K), ret_t(1.0)));

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(one_m_w, w);
  const auto ueps
      = fmin(-1.0, 2.0 * (factor + log(a) - log1p(square(neg_v) * y)) + LOG_PI);
  const auto K_small
      = (sqrt_y * sqrt(-(ueps - sqrt(-2.0 * ueps - 2.0))) - a * wdash) / a;
  const auto K_large = sqrt_y / a - wdash;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));

  // Depending on the Ks use formula for small reaction times or large
  // reaction times (see Navarro & Fuss, 2009)
  if (K_large_value > 4 * K_small_value) {
    const auto neg_vy = neg_v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; --k) {
      auto r_k_1 = a * (2.0 * k + one_m_w);
      auto x_1 = r_k_1 - neg_vy;
      auto x_over_sqrt_y_1 = x_1 / sqrt_y;
      auto d_k_1 = std_normal_lpdf(r_k_1 / sqrt_y);
      auto temp_1 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_1)));
      auto exp_d_k_1 = exp(d_k_1);
      auto ans_1 = -temp_1 * neg_vy - sqrt_y * exp_d_k_1;

      auto x_2 = r_k_1 + neg_vy;
      auto x_over_sqrt_y_2 = x_2 / sqrt_y;
      auto temp_2 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_2)));
      auto ans_2 = temp_2 * neg_vy - sqrt_y * exp_d_k_1;
      auto r_k_2 = a * (2.0 * k + 1.0 + w);
      auto d_k_2 = std_normal_lpdf(r_k_2 / sqrt_y);

      auto x_3 = r_k_2 - neg_vy;
      auto x_over_sqrt_y_3 = x_3 / sqrt_y;
      auto temp_3 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_3)));
      auto exp_d_k_2 = exp(d_k_2);
      auto ans_3 = -temp_3 * neg_vy - sqrt_y * exp_d_k_2;

      auto x_4 = r_k_2 + neg_vy;
      auto x_over_sqrt_y_4 = x_4 / sqrt_y;
      auto temp_4 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_4)));
      auto ans_4 = temp_4 * neg_vy - sqrt_y * exp_d_k_2;

      ans += (ans_1 + ans_2 + ans_3 - ans_4) * (2.0 * k + one_m_w)
             + ans_3 * one_m_w;
    }
    F_k = make_finite(exp(one_m_w_a_neg_v + 0.5 * square(neg_v) * y));
    const auto summands_small_y = ans / (y * F_k);
    return -one_m_w_neg_v * cdf + summands_small_y;
  }
  ret_t ans = 0.0;
  for (auto k = K_large_value; k > 0; --k) {
    const auto kpi = k * pi();
    const auto kpia2 = square(kpi / a);
    const auto denom = square(neg_v) + kpia2;
    auto last = (square(kpi) / pow(a, 3) * (y + 2.0 / denom)) * k / denom
                * exp(-0.5 * kpia2 * y);
    ans -= last * sin(kpi * one_m_w);
  }
  const ret_t prob
      = make_finite(exp(log_probability_distribution(a, neg_v, one_m_w)));
  const auto dav = log_probability_GradAV(a, neg_v, one_m_w);
  auto dav_neg_v = dav * neg_v;
  auto prob_deriv = fabs(neg_v) == 0
                        ? ret_t(0.0)
                        : is_inf(dav_neg_v) ? NEGATIVE_INFTY : dav_neg_v * prob;
  ans = (-2.0 / a - one_m_w_neg_v) * (cdf - prob)
        + ans * (2.0 * pi() / square(a))
              * exp(-one_m_w_a_neg_v - 0.5 * square(neg_v) * y);
  return prob_deriv + ans;
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'v' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The value of the distribution
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to v
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_v(const T_y& y, const T_a& a, const T_v& v,
                               const T_w& w, T_cdf&& cdf,
                               T_err log_err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto neg_v = -v;
  const auto one_m_w = 1.0 - w;

  const auto one_m_w_a = one_m_w * a;
  const auto one_m_w_a_neg_v = one_m_w_a * neg_v;

  const auto log_y = log(y);
  const auto factor = one_m_w_a_neg_v + square(neg_v) * y / 2.0 + log_err;

  const auto log_a = log(a);
  auto K_large_value = ret_t(1.0);
  if (neg_v != 0) {
    const auto temp = -make_finite(exp(log_a - LOG_PI - 0.5 * log_y));
    const auto log_v = log(fabs(neg_v));
    auto alphK_large = fmin(exp(factor + 0.5 * (7 * LOG_PI + log_y)
                                - 2.5 * LOG_TWO - 3 * log_a - log_v),
                            1.0);
    K_large_value
        = fmax(ceil((alphK_large == 0)
                        ? ret_t(INFTY)
                        : (alphK_large == 1) ? ret_t(NEGATIVE_INFTY)
                                             : temp * inv_Phi(alphK_large)),
               ret_t(1.0));
  }

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(one_m_w, w);
  auto K_large = fabs(neg_v) / a * y - wdash;
  const auto alphK_small = factor + 0.5 * (LOG_TWO - log_y + LOG_PI);
  const auto K_small
      = (alphK_small < 0) ? sqrt_y * sqrt(-2.0 * alphK_small) / a - wdash : 0;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));
  if (K_large_value > 4 * K_small_value) {
    const auto sqrt_y = sqrt(y);
    const auto neg_vy = neg_v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; --k) {
      auto r_k_1 = a * (2.0 * k + one_m_w);
      auto d_k_1 = std_normal_lpdf(r_k_1 / sqrt_y);
      auto x_1 = r_k_1 - neg_vy;
      auto x_over_sqrt_y_1 = x_1 / sqrt_y;
      auto ans_1 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_1)));

      auto x_2 = r_k_1 + neg_vy;
      auto x_over_sqrt_y_2 = x_2 / sqrt_y;
      auto ans_2 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_2)));
      auto r_k_2 = a * (2.0 * k + 1.0 + w);
      auto d_k_2 = std_normal_lpdf(r_k_2 / sqrt_y);

      auto x_3 = r_k_2 - neg_vy;
      auto x_over_sqrt_y_3 = x_3 / sqrt_y;
      auto ans_3 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_3)));

      auto x_4 = r_k_2 + neg_vy;
      auto x_over_sqrt_y_4 = x_4 / sqrt_y;
      auto ans_4 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_4)));
      ans += -ans_1 * x_1 + ans_2 * x_2 + ans_3 * x_3 - ans_4 * x_4;
    }
    F_k = make_finite(exp(one_m_w_a_neg_v + 0.5 * square(neg_v) * y));
    const auto summands_small_y = ans / F_k;
    return (one_m_w_a + neg_vy) * cdf - summands_small_y;
  }
  ret_t ans = 0.0;
  for (auto k = K_large_value; k > 0; --k) {
    const auto kpi = k * pi();
    const auto kpia2 = square(kpi / a);
    const auto ekpia2y = exp(-0.5 * kpia2 * y);
    const auto denom = square(neg_v) + kpia2;
    const auto denomk = k / denom;
    auto last = denomk * ekpia2y / denom;
    ans -= last * sin(kpi * one_m_w);
  }
  const ret_t prob
      = make_finite(exp(log_probability_distribution(a, neg_v, one_m_w)));
  const auto dav = log_probability_GradAV(a, neg_v, one_m_w);
  auto dav_a = dav * a;
  auto prob_deriv = is_inf(dav_a) ? ret_t(NEGATIVE_INFTY) : dav_a * prob;
  ans = (-one_m_w_a + v * y) * (cdf - prob)
        + ans * 4.0 * v * pi() / square(a)
              * exp(-one_m_w_a_neg_v - 0.5 * square(neg_v) * y);
  return -(prob_deriv + ans);
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'w' (natural-scale)
 *
 * @param y The reaction time in seconds
 * @param a The boundary separation
 * @param v The relative starting point
 * @param w The drift rate
 * @param cdf The value of the distribution
 * @param log_err The log error tolerance in the computation of the number
 * of terms for the infinite sums
 * @return Gradient with respect to w
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_w(const T_y& y, const T_a& a, const T_v& v,
                               const T_w& w, T_cdf&& cdf,
                               T_err log_err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto neg_v = -v;
  const auto one_m_w = 1 - w;

  const auto one_m_w_a_neg_v = one_m_w * a * neg_v;

  const auto factor = one_m_w_a_neg_v + square(neg_v) * y / 2.0 + log_err;

  const auto log_y = log(y);
  const auto log_a = log(a);
  const auto temp = -make_finite(exp(log_a - LOG_PI - 0.5 * log_y));
  auto alphK_large
      = fmin(exp(factor + 0.5 * (LOG_PI + log_y) - 1.5 * LOG_TWO - log_a), 1.0);
  alphK_large = fmax(0.0, alphK_large);
  const auto K_large_value
      = fmax(ceil((alphK_large == 0)
                      ? ret_t(INFTY)
                      : (alphK_large == 1) ? ret_t(NEGATIVE_INFTY)
                                           : temp * inv_Phi(alphK_large)),
             ret_t(1.0));

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(one_m_w, w);
  const auto K_large = fabs(neg_v) / a * y - wdash;
  const auto lv = log1p(square(neg_v) * y);
  const auto alphK_small = factor - LOG_TWO - lv;
  const auto arg = fmin(exp(alphK_small), 1.0);
  const auto K_small
      = (arg == 0)
            ? INFTY
            : (arg == 1) ? NEGATIVE_INFTY : -sqrt_y / a * inv_Phi(arg) - wdash;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));

  if (K_large_value > 4 * K_small_value) {
    const auto sqrt_y = sqrt(y);
    const auto neg_vy = neg_v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; --k) {
      auto r_k_1 = a * (2.0 * k + one_m_w);
      auto d_k_1 = std_normal_lpdf(r_k_1 / sqrt_y);
      auto x_1 = r_k_1 - neg_vy;
      auto x_over_sqrt_y_1 = x_1 / sqrt_y;
      auto temp_1 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_1)));
      auto exp_d_k_1 = exp(d_k_1);
      auto ans_1 = -temp_1 * neg_vy - sqrt_y * exp_d_k_1;

      auto x_2 = r_k_1 + neg_vy;
      auto x_over_sqrt_y_2 = x_2 / sqrt_y;
      auto temp_2 = make_finite(exp(d_k_1 + logMill(x_over_sqrt_y_2)));
      auto ans_2 = temp_2 * neg_vy - sqrt_y * exp_d_k_1;
      auto r_k_2 = a * (2.0 * k + 1.0 + w);
      auto d_k_2 = std_normal_lpdf(r_k_2 / sqrt_y);

      auto x_3 = r_k_2 - neg_vy;
      auto x_over_sqrt_y_3 = x_3 / sqrt_y;
      auto temp_3 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_3)));
      auto exp_d_k_2 = exp(d_k_2);
      auto ans_3 = -temp_3 * neg_vy - sqrt_y * exp_d_k_2;

      auto x_4 = r_k_2 + neg_vy;
      auto x_over_sqrt_y_4 = x_4 / sqrt_y;
      auto temp_4 = make_finite(exp(d_k_2 + logMill(x_over_sqrt_y_4)));
      auto ans_4 = temp_4 * neg_vy - sqrt_y * exp_d_k_2;

      ans += (ans_1 + ans_2 + ans_3 - ans_4) * a;
    }
    F_k = make_finite(exp(one_m_w_a_neg_v + 0.5 * square(neg_v) * y));
    const auto summands_small_y = ans / (y * F_k);
    return neg_v * a * cdf - summands_small_y;
  }
  ret_t ans = 0.0;
  for (auto k = K_large_value; k > 0; --k) {
    const auto kpi = k * pi();
    const auto kpia2 = square(kpi / a);
    const auto ekpia2y = exp(-0.5 * kpia2 * y);
    const auto denom = square(neg_v) + kpia2;
    const auto denomk = k / denom;
    auto last = kpi;
    last *= denomk * ekpia2y;
    ans -= last * cos(kpi * one_m_w);
  }
  const auto evaw = exp(-one_m_w_a_neg_v - 0.5 * square(neg_v) * y);
  const ret_t prob
      = make_finite(exp(log_probability_distribution(a, neg_v, one_m_w)));

  // Calculate the probability term 'P' on log scale
  auto dav = ret_t(-1 / w);
  if (neg_v != 0) {
    auto nearly_one = ret_t(1.0 - 1.0e-6);
    const auto sign_v = (neg_v < 0) ? 1 : -1;
    const auto sign_two_va_one_minus_w = sign_v * (2.0 * neg_v * a * w);
    const auto exp_arg = exp(sign_two_va_one_minus_w);
    if (exp_arg >= nearly_one) {
      dav = -1.0 / w;
    } else {
      auto prob = LOG_TWO + log(fabs(neg_v)) + log(a) - log1m(exp_arg);
      if (neg_v < 0) {
        prob += sign_two_va_one_minus_w;
      }
      dav = -exp(prob);
    }
  }

  const auto pia2 = 2.0 * pi() / square(a);
  auto prob_deriv = dav * prob;
  ans = v * a * (cdf - prob) + ans * pia2 * evaw;
  return -(prob_deriv + ans);
}
}  // namespace internal

/**
 * Log-CDF function for the 4-parameter Wiener distribution.
 * See 'wiener_lpdf' for more comprehensive documentation.
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
inline auto wiener_lcdf_unnorm(const T_y& y, const T_a& a, const T_t0& t0,
                               const T_w& w, const T_v& v,
                               const double& precision_derivatives = 1e-4) {
  using T_partials_return = partials_return_t<T_y, T_a, T_t0, T_w, T_v>;
  using T_y_ref = ref_type_t<T_y>;
  using T_a_ref = ref_type_t<T_a>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_w_ref = ref_type_t<T_w>;
  using T_v_ref = ref_type_t<T_v>;
  using internal::GradientCalc;
  using ret_t = return_type_t<T_y, T_a, T_t0, T_w, T_v>;

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

  if (!include_summand<propto, T_y, T_a, T_t0, T_w, T_v>::value) {
    return ret_t(0.0);
  }

  static constexpr const char* function_name = "wiener4_lcdf";
  if (size_zero(y, a, t0, w, v)) {
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
  T_partials_return lcdf = 0.0;
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

    const T_partials_return log_cdf
        = internal::estimate_with_err_check<4, 0, GradientCalc::OFF,
                                            GradientCalc::OFF>(
            [](auto&&... args) {
              return internal::wiener4_distribution<GradientCalc::OFF>(args...);
            },
            log_error_cdf - LOG_TWO, y_value - t0_value, a_value, v_value,
            w_value, log_error_absolute);

    const T_partials_return cdf = exp(log_cdf);

    lcdf += log_cdf;

    const auto new_est_err = log_cdf + log_error_derivative - LOG_FOUR;

    if (!is_constant_all<T_y>::value || !is_constant_all<T_t0>::value) {
      const auto deriv_y
          = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                              GradientCalc::ON>(
              [](auto&&... args) {
                return internal::wiener5_density<GradientCalc::ON>(args...);
              },
              new_est_err, y_value - t0_value, a_value, v_value, w_value, 0,
              log_error_absolute);

      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials)[i] = deriv_y / cdf;
      }
      if (!is_constant_all<T_t0>::value) {
        partials<2>(ops_partials)[i] = -deriv_y / cdf;
      }
    }
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                              GradientCalc::ON>(
                [](auto&&... args) {
                  return internal::wiener4_cdf_grad_a(args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / cdf;
    }
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                              GradientCalc::ON>(
                [](auto&&... args) {
                  return internal::wiener4_cdf_grad_w(args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / cdf;
    }
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener4_cdf_grad_v(y_value - t0_value, a_value, v_value,
                                         w_value, cdf, log_error_absolute)
            / cdf;
    }
  }  // for loop
  return ops_partials.build(lcdf);
}
}  // namespace math
}  // namespace stan
#endif
