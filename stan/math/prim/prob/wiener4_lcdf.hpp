#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCDF_HPP

#include <stan/math/prim/prob/wiener5_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {


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
    return ret_t(log1p(-w));
  }
  auto minus_two_va_one_minus_w = -2.0 * v * a * (1.0 - w);
  // This split prevents abort errors
  if (minus_two_va_one_minus_w < 0) {
    const auto exp_arg = exp(minus_two_va_one_minus_w);
    auto two_vaw = 2 * v * a * w;
    if (two_vaw > minus_two_va_one_minus_w) {
      return log1p(-exp_arg) - log_diff_exp(two_vaw, minus_two_va_one_minus_w);
    } else if (two_vaw < minus_two_va_one_minus_w) {
      return log1p(-exp_arg) - log_diff_exp(minus_two_va_one_minus_w, two_vaw);
    } else {
      return log1p(-exp_arg) - NEGATIVE_INFTY;
    }
  } else {
    return log1p(-exp(-minus_two_va_one_minus_w)) - log1p(-exp(2 * v * a));
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
inline auto log_probability_GradAV(const T_a& a, const T_v& v,
                                   const T_w& w) {
  using ret_t = return_type_t<T_a, T_w, T_v>;
  if (fabs(v) == 0.0) {
    return ret_t(-w);
  }
  auto nearly_one = ret_t(1.0 - 1.1 * 1.0e-5);
  //  auto nearly_one = ret_t(1.0 - std::numeric_limits<ret_t>::min());
  ret_t prob;
  if (v < 0) {
    const auto two_va_one_minus_w = (2.0 * v * a * (1.0 - w));
    const auto two_avw = 2 * a * v * w;
    const auto two_av = 2 * a * v;
    const auto exp_two_va_one_minus_w = exp(two_va_one_minus_w);
    const auto exp_two_avw = exp(two_avw);
    const auto exp_two_av = exp(two_av);
    if (((exp_two_va_one_minus_w >= nearly_one) || (exp_two_avw >= nearly_one))
        || (exp_two_av >= nearly_one)) {
      return ret_t(-w);
    }
    prob = LOG_TWO + two_va_one_minus_w - log1p(-exp_two_va_one_minus_w);
    auto log_quotient = log1p(-exp_two_avw) - log1p(-exp_two_av);
    if (log(w) > log_quotient) {
      prob += log_diff_exp(log(w), log_quotient);
      return exp(prob);
    } else {
      prob += log_diff_exp(log_quotient, log(w));
      return -exp(prob);
    }
  } else {
    const auto minus_two_va_one_minus_w = (-2.0 * v * a * (1.0 - w));
    const auto exp_minus_two_va_one_minus_w = exp(minus_two_va_one_minus_w);
    const auto minus_two_av = (-2 * a * v);
    const auto exp_minus_two_av = exp(minus_two_av);
    if ((exp_minus_two_va_one_minus_w >= nearly_one)
        || (exp_minus_two_av >= nearly_one)) {
      return ret_t(-w);
    }
    prob = LOG_TWO - log1p(-exp_minus_two_va_one_minus_w);
    ret_t log_quotient;
    if (minus_two_va_one_minus_w > minus_two_av) {
      log_quotient = log_diff_exp(minus_two_va_one_minus_w, minus_two_av)
                     - log1p(-exp_minus_two_av);
    } else if (minus_two_va_one_minus_w < minus_two_av) {
      log_quotient = log_diff_exp(minus_two_av, minus_two_va_one_minus_w)
                     - log1p(-exp_minus_two_av);
    } else {
      log_quotient = NEGATIVE_INFTY - log1p(-exp_minus_two_av);
    }
    if (log(w) > log_quotient) {
      prob += log_diff_exp(log(w), log_quotient);
      return -exp(prob);
    } else {
      prob += log_diff_exp(log_quotient, log(w));
      return exp(prob);
    }
  }
}

/**
 * Helper function. Log of Mill's ratio for the normal distribution
 *
 * @param x A scalar
 * @return Log of Mill's ratio
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
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param vn The relative starting point
 * @param wn The drift rate
 * @param wildcard This parameter space is needed for a functor. Could be
 * deleted when another solution is found
 * @param err The log error tolerance
 * @return distribution
 */
template <bool NaturalScale = false, typename T_y, typename T_a, typename T_w,
          typename T_v, typename T_wildcard, typename T_err>
inline auto wiener4_distribution(const T_y& y, const T_a& a, const T_v& vn,
                                 const T_w& wn, T_wildcard&& wildcard = 0.0,
                                 T_err&& err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto v = -vn;
  const auto w = 1 - wn;

  const auto K1 = 0.5 * (fabs(v) / a * y - w);
  const auto arg
      = fmax(0, fmin(1, exp(v * a * w + square(v) * y / 2 + err) / 2));
  const auto K2 = (arg == 0) ? INFTY
                             : (arg == 1) ? NEGATIVE_INFTY
                                          : -sqrt(y) / 2 / a * inv_Phi(arg);
  const auto K_small_value = ceil(fmax(K1, K1 + K2));

  const auto api = a / pi();
  const auto v_square = square(v);
  const auto sqrtL1 = sqrt(1 / y) * api;
  const auto sqrtL2 = sqrt(
      fmax(1.0, -2 / y * square(api)
                    * (err + log(pi() * y / 2 * (v_square + square(pi() / a)))
                       + v * a * w + v_square * y / 2)));
  const auto K_large_value = ceil(fmax(sqrtL1, sqrtL2));

  auto lg = LOG_TWO + LOG_PI - 2.0 * log(a);

  if (3 * K_small_value < K_large_value) {
    const auto sqrt_y = sqrt(y);
    const auto vy = v * y;
    auto ans = ret_t(0.0);
    ret_t fplus = NEGATIVE_INFTY;
    ret_t fminus = NEGATIVE_INFTY;
    for (auto k = K_small_value; k >= 0; k--) {
      auto rj = a * (2 * k + w);
      auto dj = std_normal_lpdf(rj / sqrt_y);
      auto pos1 = dj + logMill((rj - vy) / sqrt_y);
      auto pos2 = dj + logMill((rj + vy) / sqrt_y);
      fplus = log_sum_exp(fplus, log_sum_exp(pos1, pos2));
      rj = a * (2 * k + 2 - w);
      dj = std_normal_lpdf(rj / sqrt_y);
      auto neg1 = dj + logMill((rj - vy) / sqrt_y);
      auto neg2 = dj + logMill((rj + vy) / sqrt_y);
      fminus = log_sum_exp(fminus, log_sum_exp(neg1, neg2));
    }
    if (fplus > fminus) {
      ans = log_diff_exp(fplus, fminus);
    } else if (fplus < fminus) {
      ans = log_diff_exp(fminus, fplus);
    } else {
      ans = NEGATIVE_INFTY;
    }
    ret_t log_distribution = ans + -v * a * w - square(v) * y / 2;
    return NaturalScale ? exp(log_distribution) : log_distribution;
  } else {
    const auto log_a = log(a);
    const auto log_v = log(fabs(v));
    ret_t fplus = NEGATIVE_INFTY;
    ret_t fminus = NEGATIVE_INFTY;
    for (auto k = K_large_value; k >= 1; k--) {
      auto log_k = log(k * 1.0);
      auto k_pi = k * pi();
      auto sin_k_pi_w = sin(k_pi * w);
      if (sin_k_pi_w > 0) {
        fplus = log_sum_exp(
            fplus, log_k - log_sum_exp(2 * log_v, 2 * (log_k + LOG_PI - log_a))
                       - 0.5 * square(k_pi / a) * y + log(sin_k_pi_w));
      } else if (sin_k_pi_w < 0) {
        fminus = log_sum_exp(
            fminus, log_k - log_sum_exp(2 * log_v, 2 * (log_k + LOG_PI - log_a))
                        - 0.5 * square(k_pi / a) * y + log(-sin_k_pi_w));
      }
    }
    ret_t ans = NEGATIVE_INFTY;
    if (fplus > fminus) {
      ans = log_diff_exp(fplus, fminus);
    } else if (fplus < fminus) {
      ans = log_diff_exp(fminus, fplus);
    }
    auto summand_1 = log_probability_distribution(a, v, w);
    auto summand_2 = lg + (ans - v * a * w - 0.5 * square(v) * y);
    ret_t log_distribution = NEGATIVE_INFTY;
    if (summand_1 > summand_2) {
      log_distribution = log_diff_exp(summand_1, summand_2);
    } else if (summand_1 < summand_2) {
      log_distribution = log_diff_exp(summand_2, summand_1);
    }
    return NaturalScale ? exp(log_distribution) : log_distribution;
  }
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'a' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param vn The relative starting point
 * @param wn The drift rate
 * @param cdf The value of the distribution
 * @param err The log error tolerance
 * @return Gradient w.r.t. a
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_a(const T_y& y, const T_a& a, const T_v& vn,
                               const T_w& wn, T_cdf&& cdf,
                               T_err&& err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto v = -vn;
  const auto w = 1 - wn;

  const auto log_y = log(y);
  const auto log_a = log(a);
  auto C1
      = ret_t(LOG_TWO - log_sum_exp(2 * log(fabs(v)), 2 * (LOG_PI - log_a)));
  C1 = log_sum_exp(C1, log_y);
  const auto factor = v * a * w + square(v) * y / 2 + err;
  const auto alphK = fmin(factor + LOG_PI + log_y + log_a - LOG_TWO - C1, 0.0);
  const auto K = a / pi() / sqrt(y);
  const auto K_large_value
      = ceil(fmax(fmax(sqrt(-2 * alphK / y) * a / pi(), K), ret_t(1.0)));

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(w, 1.0 - w);
  const auto ueps = fmin(-1, 2 * (factor + log(a) - log1p(square(v) * y)) + LOG_PI);
  const auto K_small
      = (sqrt_y * sqrt(-(ueps - sqrt(-2 * ueps - 2))) - a * wdash) / a;
  const auto K_large = sqrt_y / a - wdash;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));

  if (K_large_value > 4 * K_small_value) {
    const auto vy = v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; k--) {
      auto r_k = 2 * k * a + a * w;
      auto x = r_k - vy;
      auto xsqrt_y = x / sqrt_y;
      auto d_k = std_normal_lpdf(r_k / sqrt_y);
      auto temp = fmin(exp(d_k + logMill(xsqrt_y)),
                       std::numeric_limits<ret_t>::max());
      auto temp2 = exp(d_k);
      auto temp3 = temp * (-vy) - sqrt_y * temp2;
      const auto factor = (2 * k + w);
      const auto t1 = temp3 * factor;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp3 = temp * vy - sqrt_y * temp2;
      const auto t2 = temp3 * factor;
      r_k = (2 * k + 1) * a + a * (1 - w);
      d_k = std_normal_lpdf(r_k / sqrt_y);
      x = r_k - vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp2 = exp(d_k);
      temp3 = temp * (-vy) - sqrt_y * temp2;
      const auto factor_2 = (2 * k + 2.0 - w);
      const auto t3 = -temp3 * factor_2;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp3 = temp * vy - sqrt_y * temp2;
      const auto t4 = -temp3 * factor;
      ans += (t1 + t2 + t3 + t4);
    }
    F_k = fmin(exp(v * a * w + 0.5 * square(v) * y),
               std::numeric_limits<ret_t>::max());
    const auto summands_small_y = ans / y / F_k;
    return -v * w * cdf + summands_small_y;
  } else {
    ret_t ans = NEGATIVE_INFTY;
    ans = 0.0;
    for (auto k = K_large_value; k >= 1; k--) {
      const auto kpi = k * pi();
      const auto kpia2 = square(kpi / a);
      const auto denom = 1.0 / (square(v) + kpia2);
      auto last = (square(kpi) / pow(a, 3) * (y + 2.0 * denom)) * k * denom * exp(-0.5 * kpia2 * y);
      ans += -last * sin(kpi * w);
    }
    const ret_t prob = fmin(exp(log_probability_distribution(a, v, w)),
                            std::numeric_limits<ret_t>::max());
    const auto dav = log_probability_GradAV(a, v, w);
    auto prob_deriv
        = ((fabs(v) == 0) ? ret_t(0.0)
                          : is_inf(dav * v) ? NEGATIVE_INFTY : dav * v) * prob;
    ans = (-2 / a - v * w) * (cdf - prob) + ans * (2 * pi() / square(a)) * exp(-v * a * w - 0.5 * square(v) * y);
    return prob_deriv + ans;
  }
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'v' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param vn The relative starting point
 * @param wn The drift rate
 * @param cdf The value of the distribution
 * @param err The log error tolerance
 * @return Gradient w.r.t. v
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_v(const T_y& y, const T_a& a, const T_v& vn,
                               const T_w& wn, T_cdf&& cdf,
                               T_err&& err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto v = -vn;
  const auto w = 1 - wn;
  const auto log_y = log(y);
  const auto factor = v * a * w + square(v) * y / 2 + err;

  const auto log_a = log(a);
  auto K_large_value = ret_t(1.0);
  if (v != 0) {
    const auto temp = -fmin(exp(log_a - LOG_PI - 0.5 * log_y),
                            std::numeric_limits<ret_t>::max());
    const auto log_v = log(fabs(v));
    auto alphK_large = fmin(exp(factor + 0.5 * (7 * LOG_PI + log_y)
                                - 2.5 * LOG_TWO - 3 * log_a - log_v),
                            std::numeric_limits<ret_t>::max());
    alphK_large = fmax(0.0, fmin(1.0, alphK_large));
    K_large_value
        = fmax(ceil((alphK_large == 0)
                        ? ret_t(INFTY)
                        : (alphK_large == 1) ? ret_t(NEGATIVE_INFTY)
                                             : temp * inv_Phi(alphK_large)),
               ret_t(1.0));
  }

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(w, 1.0 - w);
  auto K_large = fabs(v) / a * y - wdash;
  const auto alphK_small = factor + 0.5 * (LOG_TWO - log_y + LOG_PI);
  const auto K_small
      = (alphK_small < 0) ? sqrt_y * sqrt(-2 * alphK_small) / a - wdash : 0;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));
  if (K_large_value > 4 * K_small_value) {
    const auto sqrt_y = sqrt(y);
    const auto vy = v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; k--) {
      auto r_k = 2 * k * a + a * w;
      auto d_k = std_normal_lpdf(r_k / sqrt_y);
      auto x = r_k - vy;
      auto xsqrt_y = x / sqrt_y;
      auto temp = fmin(exp(d_k + logMill(xsqrt_y)),
                       std::numeric_limits<ret_t>::max());
      const auto factor = 2 * k + w;
      const auto factor_2 = 2 * k + 2.0 - w;
      const auto t1 = -temp * x;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      const auto t2 = temp * x;
      r_k = (2 * k + 1) * a + a * (1 - w);
      d_k = std_normal_lpdf(r_k / sqrt_y);
      x = r_k - vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      const auto t3 = temp * x;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      const auto t4 = -temp * x;
      ans += (t1 + t2 + t3 + t4);
    }
    F_k = fmin(exp(v * a * w + 0.5 * square(v) * y),
               std::numeric_limits<ret_t>::max());
    const auto summands_small_y = ans / F_k;
    return -1 * ((-w * a - v * y) * cdf + summands_small_y);
  } else {
    ret_t ans = NEGATIVE_INFTY;
    ans = 0.0;
    for (auto k = K_large_value; k >= 1; k--) {
      const auto kpi = k * pi();
      const auto factor = sin(kpi * w);
      const auto kpia2 = square(kpi / a);
      const auto ekpia2y = exp(-0.5 * kpia2 * y);
      const auto denom = 1.0 / (square(v) + kpia2);
      const auto denomk = k * denom;
      auto last = denom * denomk * ekpia2y;
      ans += -last * factor;
    }
    const ret_t prob = fmin(exp(log_probability_distribution(a, v, w)),
                            std::numeric_limits<ret_t>::max());
    const auto dav = log_probability_GradAV(a, v, w);
    auto prob_deriv = is_inf(dav * a) ? ret_t(NEGATIVE_INFTY) : dav * a;
    prob_deriv *= prob;
    ans = (-w * a - v * y) * (cdf - prob) + ans * (-2 * v) * (2 * pi() / square(a)) * exp(-v * a * w - 0.5 * square(v) * y);
    return -1 * (prob_deriv + ans);
  }
}

/**
 * Calculate derivative of the wiener4 distribution w.r.t. 'w' (natural-scale)
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param vn The relative starting point
 * @param wn The drift rate
 * @param cdf The value of the distribution
 * @param err The log error tolerance
 * @return Gradient w.r.t. w
 */
template <typename T_y, typename T_a, typename T_v, typename T_w,
          typename T_cdf, typename T_err>
inline auto wiener4_cdf_grad_w(const T_y& y, const T_a& a, const T_v& vn,
                               const T_w& wn, T_cdf&& cdf,
                               T_err&& err = log(1e-12)) {
  using ret_t = return_type_t<T_y, T_a, T_w, T_v>;
  const auto v = -vn;
  const auto w = 1 - wn;
  const auto factor = v * a * w + square(v) * y / 2 + err;

  const auto log_y = log(y);
  const auto log_a = log(a);
  const auto temp = -fmin(exp(log_a - LOG_PI - 0.5 * log_y),
                          std::numeric_limits<ret_t>::max());
  auto alphK_large
      = fmin(exp(factor + 0.5 * (LOG_PI + log_y) - 1.5 * LOG_TWO - log_a),
             std::numeric_limits<ret_t>::max());
  alphK_large = fmax(0.0, fmin(1.0, alphK_large));
  const auto K_large_value
      = fmax(ceil((alphK_large == 0)
                      ? ret_t(INFTY)
                      : (alphK_large == 1) ? ret_t(NEGATIVE_INFTY)
                                           : temp * inv_Phi(alphK_large)),
             ret_t(1.0));

  const auto sqrt_y = sqrt(y);
  const auto wdash = fmin(w, 1.0 - w);
  const auto K_large = fabs(v) / a * y - wdash;
  const auto lv = log1p(square(v) * y);
  const auto alphK_small = factor - LOG_TWO - lv;
  const auto arg
      = fmin(fmin(exp(alphK_small), std::numeric_limits<ret_t>::max()), 1.0);
  const auto K_small
      = (arg == 0)
            ? INFTY
            : (arg == 1) ? NEGATIVE_INFTY : -sqrt_y / a * inv_Phi(arg) - wdash;
  const auto K_small_value = ceil(fmax(fmax(K_small, K_large), ret_t(1.0)));

  if (K_large_value > 4 * K_small_value) {
    const auto sqrt_y = sqrt(y);
    const auto vy = v * y;
    auto ans = ret_t(0.0);
    auto F_k = ret_t(0.0);
    for (auto k = K_small_value; k >= 0; k--) {
      auto r_k = 2 * k * a + a * w;
      auto d_k = std_normal_lpdf(r_k / sqrt_y);
      auto x = r_k - vy;
      auto xsqrt_y = x / sqrt_y;
      auto temp = fmin(exp(d_k + logMill(xsqrt_y)),
                       std::numeric_limits<ret_t>::max());
      const auto factor = a;
      const auto factor_2 = -a;
      auto temp2 = exp(d_k);
      auto temp3 = temp * (-vy) - sqrt_y * temp2;
      const auto t1 = temp3 * factor;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp3 = temp * vy - sqrt_y * temp2;
      const auto t2 = temp3 * factor;
      r_k = (2 * k + 1) * a + a * (1 - w);
      d_k = std_normal_lpdf(r_k / sqrt_y);
      x = r_k - vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp2 = exp(d_k);
      temp3 = temp * (-vy) - sqrt_y * temp2;
      const auto t3 = -temp3 * factor_2;
      x = r_k + vy;
      xsqrt_y = x / sqrt_y;
      temp = fmin(exp(d_k + logMill(xsqrt_y)),
                  std::numeric_limits<ret_t>::max());
      temp3 = temp * vy - sqrt_y * temp2;
      const auto t4 = -temp3 * factor;
      ans += (t1 + t2 + t3 + t4);
    }
    F_k = fmin(exp(v * a * w + 0.5 * square(v) * y),
               std::numeric_limits<ret_t>::max());
    const auto summands_small_y = ans / y / F_k;
    return -1 * (-v * a * cdf + summands_small_y);
  } else {
    ret_t ans = NEGATIVE_INFTY;
    ans = 0.0;
    for (auto k = K_large_value; k >= 1; k--) {
      const auto kpi = k * pi();
      const auto factor = cos(kpi * w);
      const auto kpia2 = square(kpi / a);
      const auto ekpia2y = exp(-0.5 * kpia2 * y);
      const auto denom = 1.0 / (square(v) + kpia2);
      const auto denomk = k * denom;
      auto last = kpi;
      last *= denomk * ekpia2y;
      ans += -last * factor;
    }
    const auto evaw = exp(-v * a * w - 0.5 * square(v) * y);
    const ret_t prob = fmin(exp(log_probability_distribution(a, v, w)),
                            std::numeric_limits<ret_t>::max());

    // Calculate the probability term 'P' on log scale
    auto dav = ret_t(-1 / (1.0 - w));
    if (v != 0) {
      auto nearly_one = ret_t(1.0 - 1.0e-6);
      const auto sign_v = (v < 0) ? 1 : -1;
      const auto sign_two_va_one_minus_w = sign_v * (2.0 * v * a * (1.0 - w));
      const auto exp_arg = exp(sign_two_va_one_minus_w);
      if (exp_arg >= nearly_one) {
        dav = -1 / (1.0 - w);
      } else {
        auto prob = LOG_TWO + log(fabs(v)) + log(a) - log1p(-exp_arg);
        prob = (v < 0) ? prob + sign_two_va_one_minus_w : prob;
        dav = -exp(prob);
      }
    }

    const auto pia2 = 2 * pi() / square(a);
    auto prob_deriv = dav;
    prob_deriv *= prob;
    ans = -v * a * (cdf - prob) + ans * pia2 * evaw;
    return -1 * (prob_deriv + ans);
  }
}
}  // namespace internal

/**
 * Log-CDF function for the 4-parameter Wiener distribution.
 * See 'wiener_lpdf' for more comprehensive documentation
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
inline auto wiener_lcdf(const T_y& y, const T_a& a, const T_t0& t0,
                        const T_w& w, const T_v& v,
                        const double& precision_derivatives = 1e-4) {
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

  static constexpr const char* function_name = "wiener4_lcdf";
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
  T_partials_return lcdf = 0.0;
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
    const T_partials_return log_cdf
        = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                            GradientCalc::OFF>(
            [](auto&&... args) {
              return internal::wiener4_distribution<GradientCalc::OFF>(args...);
            },
            log_error_cdf - LOG_TWO, y_value - t0_value, a_value, v_value,
            w_value, 0.0, log_error_absolute);

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
