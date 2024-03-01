#ifndef STAN_MATH_PRIM_PROB_WIENER4_LCDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER4_LCDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {
namespace internal {
	

/**
 * Helper function. Calculate \f$\mathbb{e}^{x}\ if \f$x\ is not too large.
 *
 * @param x A scalar
 * @return \f$\mathbb{e}^{x}\
 */
template <typename T_x>
inline auto rexp(T_x&& x) noexcept {
  return (x <= 700) ? exp(x) : exp(700);
}

// DOKU!!
/**
 * Helper function. ?
 *
 * @param x A scalar
 * @return ?
 */
template <typename T_x>
inline auto lognormal(T_x&& x) noexcept {
  return -0.5 * x * x - LOG_SQRT_PI - 0.5 * LOG_TWO;
}

// DOKU!!
/**
 * Helper function. ?
 *
 * @param z A scalar
 * @return ?
 */
/** The original function lnnorm was written by ...
        Jean Marie Linhart
        StataCorp LP
        jlinhart@stata.com
        January 4, 2008
and later modified and adapted to Stan
*/
template <typename T_z>
inline auto lnnorm(T_z&& z) noexcept {
  using ret_t = return_type_t<T_z>;
  static const double LNNORM_MAX_X = 38;
  static const double LNNORM_MIN_X = -1e9;
  int lower;

  auto t = 0.0; 
  auto z_ = z;
  auto s = z_; 

  if (z_ == 0) {
    return (log(0.5));
  }
  if (z_ > LNNORM_MAX_X) {
    return ret_t(0);
  }
  if (z_ <= LNNORM_MIN_X) {
    return (-0.5 * z_ * z_);
  }
  if (z_ < 0) {
    z_ = -z_;
    lower = 1;
  } else
    lower = 0;

  const auto z2 = z_ * z_;
  const auto y = exp(-0.5 * z2) / SQRT_TWO_PI;
  auto n = y / z_;

  if (!(z_ > 2)) {
    z_ *= y;
    s = z_;
    for (n = 3; s != t; n += 2) {
      t = s;
      z_ *= (z2 / n);
      s += z_;
    }
    if (lower) {
      return (log(0.5 - s));
    }
    return (log(0.5 + s));
  }

  auto a1 = 2;
  auto a2 = 0;
  n = z2 + 3;
  auto p1 = 1;
  auto q1 = z_;
  auto p2 = (n - 1);
  auto q2 = n * z_;
  auto m = p1 / q1;
  t = p2 / q2;
  s = m;

  for (n += 4; m != t && s != t; n += 4) {
    a1 -= 8;
    a2 += a1;
    s = a2 * p1 + n * p2;
    p1 = p2;
    p2 = s;
    s = a2 * q1 + n * q2;
    q1 = q2;
    q2 = s;
    s = m;
    m = t;
    if (q2 > 1e30) {
      p1 /= 1e30;
      p2 /= 1e30;
      q1 /= 1e30;
      q2 /= 1e30;
    }
    t = p2 / q2;
  }
  t = lower ? log(t) - 0.5 * z2 - log(SQRT_TWO_PI) : log1p(-y * t);
  return ret_t(t);
}

// DOKU!!
/**
 * Helper function. ?
 *
 * @param x A scalar
 * @return ?
 */
template <typename T_x>
inline auto logMill(T_x&& x) noexcept {
  if (x > 1.0e5) {
    return -log(x);
  }
  const auto m = lnnorm(-x) - lognormal(x);
  return m;
}

/**
 * Calculate the 'P' term
 *
 * @tparam Distribution Whether the calculation is for the distribution
 * @tparam GradAV Whether the calculation is for gradient w.r.t. 'a' or
 * w.r.t. 'v'
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @return 'P' term
 */
template <GradientCalc Distribution = GradientCalc::OFF, GradientCalc GradAV = GradientCalc::OFF,
          GradientCalc GradW = GradientCalc::OFF, typename T_a, typename T_w, typename T_v>
inline auto logP(const T_a& a, const T_v& v, const T_w& w) noexcept {
  auto em1 = 1.0 - 1.0e-6;
  if (Distribution) {
    if (fabs(v) == 0.0) {
      return log1p(-w);
    }
    auto e = (-2.0 * v * a * (1.0 - w));
	auto prob = e; // to have an initializer.
    if (e < 0) {
      auto tt = exp(e);
      if (tt >= em1) {
        return log1p(-w);
      }
      auto two_vaw = 2 * v * a * w;
      if (two_vaw > e) {
        tt = log1p(-tt) - log_diff_exp(two_vaw, e);
      } else if (two_vaw < e) {
        tt = log1p(-tt) - log_diff_exp(e, two_vaw);
      } else {
        tt = log1p(-tt) - NEGATIVE_INFTY;
      }
      prob = tt;
    } else {
      auto tt = exp(-e);
      if (tt >= em1)
        return log1p(-w);
      tt = log1p(-tt) - log1p(-exp(2 * v * a));
      prob = tt;
    }
    return prob;
  }
  if (GradAV) {
    if (fabs(v) == 0.0) {
      return -w;
    }
    em1 = 1.0 - 1.1 * 1.0e-5;
	auto tt = 0.0; // to have an initializer.
    if (v < 0) {
      const auto emw = (2.0 * v * a * (1.0 - w));
      const auto ew = 2 * a * v * w;
      const auto e = 2 * a * v;
      const auto eemw = exp(emw);
      const auto eew = exp(ew);
      const auto ee = exp(e);
      if (((eemw >= em1) || (eew >= em1)) || (ee >= em1)) {
        return -w;
      }
      tt = LOG_TWO + emw - log1p(-eemw);
      auto temp = log1p(-eew) - log1p(-ee);
      auto sign_tt = sign(tt);
      if (log(w) > temp) {
        tt += log_diff_exp(log(w), temp);
        tt = exp(tt);
      } else {
        tt += log_diff_exp(temp, log(w));
        tt = -exp(tt);
      }
    } else {
      const auto emw = (-2.0 * v * a * (1.0 - w)), e = (-2 * a * v);
      const auto eemw = exp(emw), ee = exp(e);
      if ((eemw >= em1) || (ee >= em1)) {
        return -w;
      }
      tt = LOG_TWO - log1p(-eemw);
      auto temp = e; // to have an initializer.
      if (emw > e) {
        temp = log_diff_exp(emw, e) - log1p(-ee);
      } else if (emw < e) {
        temp = log_diff_exp(e, emw) - log1p(-ee);
      } else {
        temp = NEGATIVE_INFTY - log1p(-ee);
      }
      if (log(w) > temp) {
        tt += log_diff_exp(log(w), temp);
        tt = -exp(tt);
      } else {
        tt += log_diff_exp(temp, log(w));
        tt = exp(tt);
      }
    }
    return (is_inf(tt)) ? NEGATIVE_INFTY : tt;
  }
  if (GradW) {
    if (fabs(v) == 0.0) {
      return -1 / (1.0 - w);
    }
    const auto sign_e = (v < 0) ? 1 : -1;
    const auto e = sign_e * (2.0 * v * a * (1.0 - w));
    const auto ee = exp(e);
    if (ee >= em1) {
      return -1 / (1.0 - w);
    }
    auto temp = LOG_TWO + log(fabs(v)) + log(a) - log1p(-ee);
    temp = (v < 0) ? temp + e : temp;
    return -exp(temp);
  }
}

/**
 * Calculate 'K_small', the number of terms needed for small y
 *
 * @tparam Distribution Whether the calculation is for the distribution
 * @tparam GradA Whether the calculation is for gradient w.r.t. 'a'
 * @tparam GradV Whether the calculation is for gradient w.r.t. 'v'
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param precision The error tolerance
 * @return 'K_small' term
 */
template <GradientCalc Distribution = GradientCalc::OFF, GradientCalc GradA = GradientCalc::OFF,
          GradientCalc GradV = GradientCalc::OFF, GradientCalc GradW = GradientCalc::OFF, typename T_y, typename T_a,
          typename T_v, typename T_w, typename T_prec>
inline auto K_Small(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                      T_prec&& precision) noexcept {
  if (Distribution) {
    const auto K1 = 0.5 * (fabs(v) / a * y - w);
    const auto arg = fmax(
        0, fmin(1, exp(v * a * w + square(v) * y / 2 + (precision)) / 2));
    const auto K2 = (arg == 0) ? INFTY
                                 : (arg == 1) ? NEGATIVE_INFTY
                                              : -sqrt(y) / 2 / a * inv_Phi(arg);
    return ceil(fmax(K1, K1 + K2));
  }
  const auto sqrt_y = sqrt(y);
  const auto factor = v * a * w + square(v) * y / 2 + precision;
  const auto wdash = fmin(w, 1.0 - w);
  auto K_small = 0.0; // to have initializer
  auto K_large = fabs(v) / a * y - wdash;
  if (GradA) {
    const auto lv = log1p(square(v) * y);
    K_large = sqrt_y / a - wdash;
    const auto ueps = fmin(-1, 2 * (factor + log(a) - lv) + LOG_PI);
    K_small = (sqrt_y * sqrt(-(ueps - sqrt(-2 * ueps - 2))) - a * wdash) / a;
  }
  if (GradV) {
    const auto log_y = log(y);
    const auto alphK = factor + 0.5 * (LOG_TWO - log_y + LOG_PI);
    K_small = (alphK < 0) ? sqrt_y * sqrt(-2 * alphK) / a - wdash : 0;
  }
  if (GradW) {
    const auto lv = log1p(square(v) * y);
    const auto alphK = factor - LOG_TWO - lv;
    const auto arg = fmin(rexp(alphK), 1.0);
    K_small = (arg == 0) ? INFTY
                         : (arg == 1) ? NEGATIVE_INFTY
                                      : -sqrt_y / a * inv_Phi(arg) - wdash;
  }
  return ceil(fmax(fmax(K_small, K_large), 1.0));
}

/**
 * Calculate 'K_large', the number of terms needed for large y
 *
 * @tparam Distribution Whether the calculation is for the distribution
 * @tparam GradA Whether the calculation is for gradient w.r.t. 'a'
 * @tparam GradV Whether the calculation is for gradient w.r.t. 'v'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param precision The error tolerance
 * @return 'K_large' term
 */
template <GradientCalc Distribution = GradientCalc::OFF, GradientCalc GradA = GradientCalc::OFF,
          GradientCalc GradV = GradientCalc::OFF, typename T_y, typename T_a, typename T_v,
          typename T_w, typename T_prec>
inline auto K_Large(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                      T_prec&& precision) noexcept {
  if (Distribution) {
    const auto api = a / pi();
    const auto v_square = square(v);
    const auto sqrtL1 = sqrt(1 / y) * api;
    const auto sqrtL2 = sqrt(fmax(
        1.0,
        -2 / y * square(api)
            * (precision + log(pi() * y / 2 * (v_square + square(pi() / a)))
               + v * a * w + v_square * y / 2)));
    return ceil(fmax(sqrtL1, sqrtL2));
  }
  const auto log_y = log(y);
  const auto log_a = log(a);
  const auto factor = v * a * w + square(v) * y / 2 + precision;
  if (GradA) {
    const auto K = a / pi() / sqrt(y);
    auto C1 = LOG_TWO - log_sum_exp(2 * log(fabs(v)), 2 * (LOG_PI - log_a));
    C1 = log_sum_exp(C1, log_y);
    const auto alphK
        = fmin(factor + LOG_PI + log_y + log_a - LOG_TWO - C1, 0.0);
    return ceil(fmax(fmax(sqrt(-2 * alphK / y) * a / pi(), K), 1.0));
  }
  if (GradV && v == 0) {
    return 1.0;
  } else {
    const auto temp = -rexp(log_a - LOG_PI - 0.5 * log_y);
    auto alphK = 0.0; //to have initializer
    if (GradV) {
      const auto log_v = log(fabs(v));
      alphK = rexp(factor + 0.5 * (7 * LOG_PI + log_y) - 2.5 * LOG_TWO
                           - 3 * log_a - log_v);
      alphK = fmax(0.0, fmin(1.0, alphK));
    } else {
      alphK = rexp(factor + 0.5 * (LOG_PI + log_y) - 1.5 * LOG_TWO
                           - log_a);
      alphK = fmax(0.0, fmin(1.0, alphK));
    }
    return fmax(ceil((alphK == 0) ? INFTY
                                  : (alphK == 1) ? NEGATIVE_INFTY
                                                 : temp * inv_Phi(alphK)),
                1.0);
  }
}

/**
 * Calculate terms of the sum for small y
 *
 * @tparam Distribution Whether the calculation is for the distribution
 * @tparam GradA Whether the calculation is for gradient w.r.t. 'a'
 * @tparam GradV Whether the calculation is for gradient w.r.t. 'v'
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param K The number of summands
 * @return terms of the sum
 */
template <GradientCalc Distribution = GradientCalc::OFF, GradientCalc GradA = GradientCalc::OFF,
          GradientCalc GradV = GradientCalc::OFF, GradientCalc GradW = GradientCalc::OFF, typename T_y, typename T_a,
          typename T_v, typename T_w, typename T_K>
inline auto logFs(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                    T_K&& K) noexcept {
  const auto sqrt_y = sqrt(y);
  const auto vy = v * y;
  auto ans = 0.0;
  if (Distribution) {
    auto fplus = NEGATIVE_INFTY;
    auto fminus = NEGATIVE_INFTY;
    for (auto k = K; k >= 0; k--) {
      auto rj = a * (2 * k + w);
      auto dj = lognormal(rj / sqrt_y);
      auto pos1 = dj + logMill((rj - vy) / sqrt_y);
      auto pos2 = dj + logMill((rj + vy) / sqrt_y);
      fplus = log_sum_exp(fplus, log_sum_exp(pos1, pos2));
      rj = a * (2 * k + 2 - w);
      dj = lognormal(rj / sqrt_y);
      auto neg1 = dj + logMill((rj - vy) / sqrt_y);
      auto neg2 = dj + logMill((rj + vy) / sqrt_y);
      fminus = log_sum_exp(fminus, log_sum_exp(neg1, neg2));
    }
    //		auto ans_sign;
    //		std::forward_as_tuple(ans, ans_sign)
    //        = stan::math::log_sum_exp_signed(fplus, 1, fminus, -1);
    if (fplus > fminus) {
      ans = log_diff_exp(fplus, fminus);
    } else if (fplus < fminus) {
      ans = log_diff_exp(fminus, fplus);
    } else {
      ans = NEGATIVE_INFTY;
    }
    return ans + -v * a * w - square(v) * y / 2;
  }
  auto Fj = 0.0;
  for (auto k = K; k >= 0; k--) {
    auto rj = 2 * k * a + a * w;
    auto dj = lognormal(rj / sqrt_y);
    auto x = rj - vy;
    auto xsqrt_y = x / sqrt_y;
    auto temp = rexp(dj + logMill(xsqrt_y));
    auto temp2 = 0.0; //to have initializer
    auto temp3 = 0.0; //to have initializer
    auto t1 = 0.0; //to have initializer
    auto t2 = 0.0; //to have initializer
    auto t3 = 0.0; //to have initializer
    auto t4 = 0.0; //to have initializer
    const auto factor = GradW ? a : (2 * k + w);
    const auto factor_2 = GradW ? -a : (2 * k + 2.0 - w);
    if (GradA || GradW) {
      temp2 = exp(dj);
      temp3 = temp * (-vy) - sqrt_y * temp2;
      t1 = temp3 * factor;
    }
    if (GradV) {
      t1 = -temp * x;
    }
    x = rj + vy;
    xsqrt_y = x / sqrt_y;
    temp = rexp(dj + logMill(xsqrt_y));
    if (GradA || GradW) {
      temp3 = temp * vy - sqrt_y * temp2;
      t2 = temp3 * factor;
    }
    if (GradV) {
      t2 = temp * x;
    }
    rj = (2 * k + 1) * a + a * (1 - w);
    dj = lognormal(rj / sqrt_y);
    x = rj - vy;
    xsqrt_y = x / sqrt_y;
    temp = rexp(dj + logMill(xsqrt_y));
    if (GradA || GradW) {
      temp2 = exp(dj);
      temp3 = temp * (-vy) - sqrt_y * temp2;
      t3 = -temp3 * factor_2;
    }
    if (GradV) {
      t3 = temp * x;
    }
    x = rj + vy;
    xsqrt_y = x / sqrt_y;
    temp = rexp(dj + logMill(xsqrt_y));
    if (GradA || GradW) {
      temp3 = temp * vy - sqrt_y * temp2;
      t4 = -temp3 * factor;
    }
    if (GradV) {
      t4 = -temp * x;
    }
    ans += (t1 + t2 + t3 + t4);
  }
  Fj = rexp(v * a * w + 0.5 * square(v) * y);
  return (GradA || GradW) ? ans / y / Fj : ans / Fj;
}

/**
 * Calculate terms of the sum for large y
 *
 * @tparam Distribution Whether the calculation is for the distribution
 * @tparam GradA Whether the calculation is for gradient w.r.t. 'a'
 * @tparam GradV Whether the calculation is for gradient w.r.t. 'v'
 * @tparam GradW Whether the calculation is for gradient w.r.t. 'w'
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param K The number of summands
 * @return terms of the sum
 */
template <GradientCalc Distribution = GradientCalc::OFF, GradientCalc GradA = GradientCalc::OFF,
          GradientCalc GradV = GradientCalc::OFF, GradientCalc GradW = GradientCalc::OFF, typename T_y, typename T_a,
          typename T_v, typename T_w, typename T_K, typename T_cdf>
inline auto logFl(const T_y& y, const T_a& a, const T_v& v, const T_w& w,
                    T_K&& K, T_cdf&& cdf = 0.0) noexcept {
  auto ans = NEGATIVE_INFTY;
  if (Distribution) {
    const auto log_a = log(a);
    const auto log_v = log(fabs(v));
    auto fplus = NEGATIVE_INFTY;
    auto fminus = NEGATIVE_INFTY;
    //		int ans_sign = -1;
    for (auto k = K; k >= 1; k--) {
      auto temp0 = log(k * 1.0);
      auto temp1 = k * pi();
      auto temp2 = temp1 * w;
      auto check = sin(temp2);
      //			int check_sign = sign(check);
      //			auto summand_2 = temp0 - log_sum_exp(2 *
      //log_v, 2 * (temp0 + LOG_PI - log_a)) - 0.5 * square(temp1 / a) * y +
      //log(fabs(check)); 			std::forward_as_tuple(ans, ans_sign) 			=
      //stan::math::log_sum_exp_signed(ans, ans_sign, summand_2, check_sign);
      if (check > 0) {
        fplus = log_sum_exp(
            fplus, temp0 - log_sum_exp(2 * log_v, 2 * (temp0 + LOG_PI - log_a))
                       - 0.5 * square(temp1 / a) * y + log(check));
      } else if (check < 0) {
        fminus = log_sum_exp(
            fminus, temp0 - log_sum_exp(2 * log_v, 2 * (temp0 + LOG_PI - log_a))
                        - 0.5 * square(temp1 / a) * y + log(-check));
      }
    }
    if (fplus > fminus) {
      ans = log_diff_exp(fplus, fminus);
    } else if (fplus < fminus) {
      ans = log_diff_exp(fminus, fplus);
    }
    return (ans - v * a * w - 0.5 * square(v) * y);
  }
  ans = 0.0;
  for (auto k = K; k >= 1; k--) {
    const auto kpi = k * pi();
    const auto factor = (GradA || GradV) ? sin(kpi * w) : cos(kpi * w);
    const auto kpia2 = square(kpi / a);
    const auto ekpia2y = exp(-0.5 * kpia2 * y);
    const auto denom = 1.0 / (square(v) + kpia2);
    const auto denomk = k * denom;
    auto last = GradA ? square(kpi) / pow(a, 3) * (y + 2.0 * denom)
                        : GradV ? denom : kpi;
    last *= denomk * ekpia2y;
    ans += -last * factor;
  }
  const auto evaw = exp(-v * a * w - 0.5 * square(v) * y);
  const auto temp = rexp(logP<GradientCalc::ON>(a, v, w));
  const auto GradAV = (GradA || GradV)? GradientCalc::ON : GradientCalc::OFF;
  const auto dav = logP<GradientCalc::OFF, GradAV, GradW>(a, v, w);
  const auto pia2 = 2 * pi() / square(a);
  auto temp_deriv
      = GradA
            ? ((fabs(v) == 0) ? 0 : is_inf(dav * v) ? NEGATIVE_INFTY : dav * v)
            : GradV ? is_inf(dav * a) ? NEGATIVE_INFTY : dav * a : dav;
  temp_deriv *= temp;
  ans = GradA ? (-2 / a - v * w) * (cdf - temp) + ans * pia2 * evaw
              : GradV ? (-w * a - v * y) * (cdf - temp)
                            + ans * (-2 * v) * pia2 * evaw
                      : -v * a * (cdf - temp) + ans * pia2 * evaw;
  return temp_deriv + ans;
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
 * @param err The log error tolerance
 * @return distribution
 */
template <bool NaturalScale = false, typename T_y,
          typename T_a, typename T_w, typename T_v, typename T_wildcard>
inline auto wiener4_distribution(const T_y& y, const T_a& a, const T_v& vn,
                                   const T_w& wn, T_wildcard&& wildcard = 0,
                                   const double& err = log(1e-12)) noexcept {
  auto ans = 0.0; //to have initializer
  const auto v = -vn;
  const auto w = 1 - wn;

  if (is_inf(y)) {
    logP<GradientCalc::ON>(a, v, w);
  }

  auto K_small = K_Small<GradientCalc::ON>(y, a, v, w, err);
  auto K_large = K_Large<GradientCalc::ON>(y, a, v, w, err);
  auto lg = LOG_TWO + LOG_PI - 2.0 * log(a);

  if (3 * K_small < K_large) {
    ans = logFs<GradientCalc::ON>(y, a, v, w, K_small);
  } else {
    auto summand_1 = logP<GradientCalc::ON>(a, v, w);
    auto summand_2 = lg + logFl<GradientCalc::ON>(y, a, v, w, K_large, 0);
    //		Scalar ans_sign;
    //		std::forward_as_tuple(ans, ans_sign)
    //       = stan::math::log_sum_exp_signed(summand_1, 1, summand_2, -1);
    if (summand_1 > summand_2) {
      ans = log_diff_exp(summand_1, summand_2);
    } else if (summand_1 < summand_2) {
      ans = log_diff_exp(summand_2, summand_1);
    } else {
      ans = NEGATIVE_INFTY;
    }
  }
  return NaturalScale ? exp(ans) : ans;
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
template <typename T_y, typename T_a, typename T_v,
          typename T_w, typename T_cdf>
inline auto wiener4_cdf_grad_a(const T_y& y, const T_a& a, const T_v& vn,
                                 const T_w& wn, T_cdf&& cdf,
                                 const double& err = log(1e-12)) noexcept {
  const auto v = -vn;
  const auto w = 1 - wn;
  const auto K_large = K_Large<GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, err);
  const auto K_small = K_Small<GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, err);
  if (K_large > 4 * K_small) {
    return -v * w * cdf + logFs<GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, (K_small));
  } else {
    return logFl<GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, K_large, cdf);
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
template <typename T_y, typename T_a, typename T_v,
          typename T_w, typename T_cdf>
inline auto wiener4_cdf_grad_v(const T_y& y, const T_a& a, const T_v& vn,
                                 const T_w& wn, T_cdf&& cdf,
                                 const double& err = log(1e-12)) noexcept {
  const auto v = -vn;
  const auto w = 1 - wn;
  const auto K_large = K_Large<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, err);
  const auto K_small = K_Small<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, err);
  if (K_large > 4 * K_small) {
    return -1
           * ((-w * a - v * y) * cdf
              + logFs<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, (K_small)));
  } else {
    return -1 * logFl<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, (K_large), cdf);
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
template <typename T_y, typename T_a, typename T_v,
          typename T_w, typename T_cdf>
inline auto wiener4_cdf_grad_w(const T_y& y, const T_a& a, const T_v& vn,
                                 const T_w& wn, T_cdf&& cdf,
                                 const double& err = log(1e-12)) noexcept {
  const auto v = -vn;
  const auto w = 1 - wn;
  const auto K_large = K_Large(y, a, v, w, err);
  const auto K_small
      = K_Small<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, err);
  if (K_large > 4 * K_small) {
    return -1
           * (-v * a * cdf
              + logFs<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w,
                                                         K_small));
  } else {
    return -1
           * logFl<GradientCalc::OFF, GradientCalc::OFF, GradientCalc::OFF, GradientCalc::ON>(y, a, v, w, K_large,
                                                      cdf);
  }
}
}  // namespace internal

/**
 * Log-CDF function for the 4-parameter Wiener distribution.
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
 * @return The log of the Wiener first passage time distribution with
 *  the specified arguments for upper boundary responses
 */
template <bool propto = false, typename T_y, typename T_a, typename T_t0,
          typename T_w, typename T_v>
inline auto wiener4_lcdf(
    const T_y& y, const T_a& a, const T_t0& t0, const T_w& w, const T_v& v,
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
  const auto log_error_absolute = log(1e-12);
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
    const auto log_cdf
        = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                                       GradientCalc::OFF>(
            [](auto&&... args) {
              return internal::wiener4_distribution<GradientCalc::OFF>(
                  args...);
            },
            log_error_cdf - LOG_TWO, y_value - t0_value, a_value, v_value,
            w_value, 0, log_error_absolute);

    const auto cdf = exp(log_cdf);

    lcdf += log_cdf;

    const auto new_est_err
        = log_cdf + log_error_derivative - LOG_FOUR;

    const auto deriv_y
        = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                            GradientCalc::ON>(
            [](auto&&... args) {
              return internal::wiener5_density<GradientCalc::ON>(
                  args...);
            },
            new_est_err, y_value - t0_value, a_value, v_value, w_value, 0,
            log_error_absolute);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[i] = deriv_y / cdf;
    }
    if (!is_constant_all<T_a>::value) {
      partials<1>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                              GradientCalc::ON>(
                [](auto&&... args) {
                  return internal::wiener4_cdf_grad_a(
                      args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / cdf;
    }
    if (!is_constant_all<T_t0>::value) {
      partials<2>(ops_partials)[i] = -deriv_y / cdf;
    }
    if (!is_constant_all<T_w>::value) {
      partials<3>(ops_partials)[i]
          = internal::estimate_with_err_check<5, 0, GradientCalc::OFF,
                                              GradientCalc::ON>(
                [](auto&&... args) {
                  return internal::wiener4_cdf_grad_w(
                      args...);
                },
                new_est_err, y_value - t0_value, a_value, v_value, w_value, cdf,
                log_error_absolute)
            / cdf;
    }
    if (!is_constant_all<T_v>::value) {
      partials<4>(ops_partials)[i]
          = internal::wiener4_cdf_grad_v(
                y_value - t0_value, a_value, v_value, w_value, cdf)
            / cdf;
    }
  }  // for loop
  return ops_partials.build(lcdf);
}
}  // namespace math
}  // namespace stan
#endif
