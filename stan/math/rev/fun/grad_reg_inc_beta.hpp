#ifndef STAN_MATH_REV_FUN_GRAD_REG_INC_BETA_HPP
#define STAN_MATH_REV_FUN_GRAD_REG_INC_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <cmath>
#include <limits>
#include <utility>

namespace stan {
namespace math {

namespace internal {

/**
 * Hessian-vector product H v of the regularized incomplete beta function
 * I_x(a, b), H its Hessian in (a, b), by the power series
 *
 *   I_x(a, b) = N S,  N = x^a / B(a, b),  S = sum_n P_n / (a + n),
 *   P_n = prod_{k <= n} (k - b) x / k.
 *
 * Differentiating I = N S twice,
 *
 *   H v = N (D_v grad S + grad S D_v log N + grad log N D_v S
 *            + S (grad log N D_v log N + D_v grad log N)),
 *
 * and a diagonal entry is
 *
 *   H_jj = N (d2S + 2 dlog_N dS + S (dlog_N^2 + d2log_N)),
 *
 * d and d2 the first and second derivatives in shape j. The a derivatives
 * of each term t = P_n / (a + n) are closed form, and its b derivatives need
 * those of P_n, carried by recurrence.
 *
 * Names: for <code>G1</code> and <code>G2</code> only the requested shape
 * varies, and dX and d2X are the first and second derivatives of X in it. For
 * <code>BOTH</code>, X_a and X_b are the partials of X and Dv_X its
 * derivative along v.
 *
 * @tparam grads <code>G1</code> for H_aa v_a, <code>G2</code> for H_bb v_b,
 *   <code>BOTH</code> for H v
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param v_a a component of v; not used for <code>G2</code>
 * @param v_b b component of v; not used for <code>G1</code>
 * @param digamma_a digamma(a); not used for <code>G2</code>
 * @param digamma_b digamma(b); not used for <code>G1</code>
 * @param digamma_ab digamma(a + b)
 * @return H_aa v_a for <code>G1</code>, H_bb v_b for <code>G2</code>,
 *   ((H v)_a, (H v)_b) for <code>BOTH</code>
 */
template <GRAD_RIB grads>
inline auto inc_beta_hvp_power_series(double a, double b, double x, double v_a,
                                      double v_b, double digamma_a,
                                      double digamma_b, double digamma_ab) {
  constexpr int max_iterations = 600;
  const double tolerance = std::numeric_limits<double>::epsilon() / a;
  const double log_x = std::log(x);
  const double N = std::exp(a * log_x - lbeta(a, b));
  // boost's trigamma: stan::math::trigamma(double) is good to about 1e-9 below
  // 5, and psi1(a + b) enters H_ab alone, where that error does not cancel
  const double trigamma_ab = boost::math::trigamma(a + b, boost_policy_t<>());
  // n = 0: P_0 = 1 and the term is 1 / a, with a derivatives -1 / a^2 and
  // 2 / a^3, and it does not depend on b
  double product = 1.0;
  double S = 1.0 / a;
  if constexpr (grads == GRAD_RIB::G1) {
    double dS = -S / a;
    double d2S = 2.0 * S / (a * a);
    for (int n = 1; n <= max_iterations; ++n) {
      product *= (n - b) * x / n;
      // dt/da = -t / (a + n), d2t/da2 = 2 t / (a + n)^2
      const double inv_apn = 1.0 / (a + n);
      const double term = product * inv_apn;
      S += term;
      dS -= term * inv_apn;
      d2S += 2.0 * term * inv_apn * inv_apn;
      // for integer b the term and its a derivatives are 0 from n = b on
      if (std::fabs(term) <= tolerance) {
        break;
      }
    }
    const double dlog_N = log_x - digamma_a + digamma_ab;
    const double d2log_N
        = trigamma_ab - boost::math::trigamma(a, boost_policy_t<>());
    return v_a * N
           * (d2S + 2.0 * dlog_N * dS + S * (dlog_N * dlog_N + d2log_N));
  } else if constexpr (grads == GRAD_RIB::G2) {
    double dproduct = 0.0;
    double d2product = 0.0;
    double dS = 0.0;
    double d2S = 0.0;
    for (int n = 1; n <= max_iterations; ++n) {
      // P_n = f P_{n-1} with f = (n - b) x / n and df/db = -x / n
      const double x_div_n = x / n;
      const double factor = (n - b) * x_div_n;
      d2product = factor * d2product - 2.0 * x_div_n * dproduct;
      dproduct = factor * dproduct - x_div_n * product;
      product *= factor;
      const double inv_apn = 1.0 / (a + n);
      const double term = product * inv_apn;
      const double dterm = dproduct * inv_apn;
      const double d2term = d2product * inv_apn;
      S += term;
      dS += dterm;
      d2S += d2term;
      // for integer b the term is 0 from n = b on while its b derivatives are
      // not, so every carried series has to be small before stopping
      if (std::fabs(term) <= tolerance && std::fabs(dterm) <= tolerance
          && std::fabs(d2term) <= tolerance) {
        break;
      }
    }
    const double dlog_N = digamma_ab - digamma_b;
    const double d2log_N
        = trigamma_ab - boost::math::trigamma(b, boost_policy_t<>());
    return v_b * N
           * (d2S + 2.0 * dlog_N * dS + S * (dlog_N * dlog_N + d2log_N));
  } else {
    double product_b = 0.0;
    double product_bb = 0.0;
    double S_a = -S / a;
    double S_b = 0.0;
    double Dv_S_a = 2.0 * v_a * S / (a * a);
    double Dv_S_b = 0.0;
    for (int n = 1; n <= max_iterations; ++n) {
      // P_n = f P_{n-1} with f = (n - b) x / n and df/db = -x / n
      const double x_div_n = x / n;
      const double factor = (n - b) * x_div_n;
      product_bb = factor * product_bb - 2.0 * x_div_n * product_b;
      product_b = factor * product_b - x_div_n * product;
      product *= factor;
      const double inv_apn = 1.0 / (a + n);
      const double term = product * inv_apn;
      const double term_b = product_b * inv_apn;
      // d2t/da2 = 2 t / (a + n)^2, d2t/dadb = -t_b / (a + n),
      // d2t/db2 = P''_n / (a + n)
      const double term_aa = 2.0 * term * inv_apn * inv_apn;
      const double term_ab = -term_b * inv_apn;
      const double term_bb = product_bb * inv_apn;
      S += term;
      S_a -= term * inv_apn;
      S_b += term_b;
      Dv_S_a += v_a * term_aa + v_b * term_ab;
      Dv_S_b += v_a * term_ab + v_b * term_bb;
      // for integer b the term is 0 from n = b on while its b derivatives are
      // not, so every carried series has to be small before stopping
      if (std::fabs(term) <= tolerance && std::fabs(term_b) <= tolerance
          && std::fabs(term_bb) <= tolerance) {
        break;
      }
    }
    const double log_N_a = log_x - digamma_a + digamma_ab;
    const double log_N_b = digamma_ab - digamma_b;
    const double Dv_log_N = log_N_a * v_a + log_N_b * v_b;
    const double Dv_S = S_a * v_a + S_b * v_b;
    // the Hessian of log N is
    // [[psi1(a+b) - psi1(a), psi1(a+b)], [psi1(a+b), psi1(a+b) - psi1(b)]]
    const double Dv_log_N_a
        = (trigamma_ab - boost::math::trigamma(a, boost_policy_t<>())) * v_a
          + trigamma_ab * v_b;
    const double Dv_log_N_b
        = trigamma_ab * v_a
          + (trigamma_ab - boost::math::trigamma(b, boost_policy_t<>())) * v_b;
    return std::pair{N
                         * (Dv_S_a + S_a * Dv_log_N + log_N_a * Dv_S
                            + S * (log_N_a * Dv_log_N + Dv_log_N_a)),
                     N
                         * (Dv_S_b + S_b * Dv_log_N + log_N_b * Dv_S
                            + S * (log_N_b * Dv_log_N + Dv_log_N_b))};
  }
}

/**
 * Hessian-vector product H v of the regularized incomplete beta function
 * I_x(a, b), H its Hessian in (a, b), by the continued fraction
 * https://dlmf.nist.gov/8.17.E22, I_x(a, b) = Q h with
 * Q = x^a (1 - x)^b / (a B(a, b)) and h the continued fraction, evaluated
 * with the modified Lentz method. Differentiating I = Q h twice,
 *
 *   H = I (grad log I grad log I^T + Hess log Q + Hess log h).
 *
 * log h is the sum over the Lentz steps of log C and log D, so its
 * derivatives are sums of theirs. With d the partial numerator, one step
 * takes C' = 1 + d / C and D' = 1 / (1 + d D), and
 *
 *   grad log C' = (grad d - d grad log C) / (C + d),
 *   grad log D' = -(grad d + d grad log D) D D'.
 *
 * Differentiating those once more along a direction v,
 *
 *   D_v grad log C' = (Hess d v - D_v d grad log C - d D_v grad log C)
 *                     / (C + d) - grad log C' D_v log(C C'),
 *   D_v log(C C') = (C D_v log C + D_v d) / (C + d),
 *   D_v grad log D' = -(Hess d v + D_v d grad log D + d D_v grad log D) D D'
 *                     + grad log D' D_v log(D D'),
 *   D_v log(D D') = D_v log D - D D' (D_v d + d D_v log D).
 *
 * <code>G1</code> and <code>G2</code> take v = (1, 0) and (0, 1), so their
 * D_v are second derivatives in one shape.
 *
 * Names: for <code>G1</code> and <code>G2</code> only the requested shape
 * varies, and dX and d2X are the first and second derivatives of X in it. For
 * <code>BOTH</code>, X_a and X_b are the partials of X and Dv_X its
 * derivative along v.
 *
 * @tparam grads <code>G1</code> for H_aa v_a, <code>G2</code> for H_bb v_b,
 *   <code>BOTH</code> for H v
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param v_a a component of v; not used for <code>G2</code>
 * @param v_b b component of v; not used for <code>G1</code>
 * @param digamma_a digamma(a); not used for <code>G2</code>
 * @param digamma_b digamma(b); not used for <code>G1</code>
 * @param digamma_ab digamma(a + b)
 * @return H_aa v_a for <code>G1</code>, H_bb v_b for <code>G2</code>,
 *   ((H v)_a, (H v)_b) for <code>BOTH</code>
 */
template <GRAD_RIB grads>
inline auto inc_beta_hvp_continued_fraction(double a, double b, double x,
                                            double v_a, double v_b,
                                            double digamma_a, double digamma_b,
                                            double digamma_ab) {
  constexpr int max_iterations = 300;
  const double tolerance = 3.0 * std::numeric_limits<double>::epsilon();
  const double small = std::sqrt(std::numeric_limits<double>::min());
  const double log_x = std::log(x);
  const double log1m_x = log1p(-x);
  const double Q
      = std::exp(a * log_x + b * log1m_x - std::log(a) - lbeta(a, b));
  // boost's trigamma, as in the power series
  const double trigamma_ab = boost::math::trigamma(a + b, boost_policy_t<>());

  // Partial numerator d_n = num / den over den = (a + n - 1)(a + n),
  // https://dlmf.nist.gov/8.17.E23, with its a and b derivatives:
  //   n = 2m:     num = m (b - m) x,          num_b = m x
  //   n = 2m + 1: num = -(a + m)(a + b + m) x, num_a = -(2 (a + m) + b) x,
  //               num_b = -(a + m) x, num_aa = -2 x, num_ab = -x
  // n = 1 is the odd case with m = 0, and d_bb = 0 for every n.

  // first step, n = 1, with C held at 1: D = 1 / (1 + d_1)
  const double ap1 = a + 1.0;
  const double d_1 = -(a + b) * x / ap1;
  double D = 1.0 + d_1;
  if (std::fabs(D) < small) {
    D = small;
  }
  D = 1.0 / D;
  double C = 1.0;
  double h = D;
  if constexpr (grads == GRAD_RIB::G1) {
    const double dd_1 = -(x + d_1) / ap1;
    const double d2d_1 = 2.0 * (x + d_1) / (ap1 * ap1);
    double dlog_C = 0.0;
    double d2log_C = 0.0;
    double dlog_D = -dd_1 * D;
    double d2log_D = -(d2d_1 + dlog_D * dd_1) * D;
    double dlog_h = dlog_D;
    double d2log_h = d2log_D;
    double d2log_h_prev = d2log_h;
    for (int n = 2; n <= 2 * max_iterations + 1; ++n) {
      const int m = n / 2;
      const bool odd = n % 2 == 1;
      const double a_plus_n = a + n;
      const double inv_den = 1.0 / (a_plus_n * (a_plus_n - 1.0));
      const double dden = 2.0 * a_plus_n - 1.0;
      double num;
      double dnum = 0.0;
      double d2num = 0.0;
      if (odd) {
        num = -(a + m) * (a + b + m) * x;
        dnum = -(2.0 * (a + m) + b) * x;
        d2num = -2.0 * x;
      } else {
        num = m * (b - m) * x;
      }
      const double d_n = num * inv_den;
      const double dd_n = (dnum - d_n * dden) * inv_den;
      const double d2d_n = (d2num - 2.0 * dd_n * dden - 2.0 * d_n) * inv_den;

      // one step of the modified Lentz method. In each half, the second
      // derivative takes its old-step part, the first derivative updates in
      // place, then scales dlog(C C_new) or dlog(D D_new).
      // C_new = 1 + d_n / C, with C C_new = C + d_n
      double C_new = 1.0 + d_n / C;
      double C_C_new = C + d_n;
      if (std::fabs(C_new) < small) {
        C_new = small;
        C_C_new = C * small;
      }
      const double inv_C_C_new = 1.0 / C_C_new;
      const double dlog_C_C_new = (C * dlog_C + dd_n) * inv_C_C_new;
      d2log_C = (d2d_n - dd_n * dlog_C - d_n * d2log_C) * inv_C_C_new;
      dlog_C = (dd_n - d_n * dlog_C) * inv_C_C_new;
      d2log_C -= dlog_C * dlog_C_C_new;

      // D_new = 1 / (1 + d_n D)
      double D_new = 1.0 + d_n * D;
      if (std::fabs(D_new) < small) {
        D_new = small;
      }
      D_new = 1.0 / D_new;
      const double D_D_new = D * D_new;
      const double dlog_D_D_new = dlog_D - D_D_new * (dd_n + d_n * dlog_D);
      d2log_D = -(d2d_n + dd_n * dlog_D + d_n * d2log_D) * D_D_new;
      dlog_D = -(dd_n + d_n * dlog_D) * D_D_new;
      d2log_D += dlog_D * dlog_D_D_new;

      dlog_h += dlog_C + dlog_D;
      d2log_h += d2log_C + d2log_D;
      C = C_new;
      D = D_new;
      const double delta = C_new * D_new;
      h *= delta;
      if (odd) {
        if (std::fabs(delta - 1.0) < tolerance
            && std::fabs(d2log_h - d2log_h_prev)
                   <= tolerance * std::fabs(d2log_h)) {
          break;
        }
        d2log_h_prev = d2log_h;
      }
    }
    const double dlog_I = log_x - 1.0 / a + digamma_ab - digamma_a + dlog_h;
    const double d2log_Q = trigamma_ab
                           - boost::math::trigamma(a, boost_policy_t<>())
                           + 1.0 / (a * a);
    return v_a * h * Q * (dlog_I * dlog_I + d2log_Q + d2log_h);
  } else if constexpr (grads == GRAD_RIB::G2) {
    // d_n is linear in b, so d2d_n = 0
    const double dd_1 = -x / ap1;
    double dlog_C = 0.0;
    double d2log_C = 0.0;
    double dlog_D = -dd_1 * D;
    double d2log_D = -dlog_D * dd_1 * D;
    double dlog_h = dlog_D;
    double d2log_h = d2log_D;
    double d2log_h_prev = d2log_h;
    for (int n = 2; n <= 2 * max_iterations + 1; ++n) {
      const int m = n / 2;
      const bool odd = n % 2 == 1;
      const double a_plus_n = a + n;
      const double inv_den = 1.0 / (a_plus_n * (a_plus_n - 1.0));
      double dnum;
      double num;
      if (odd) {
        dnum = -(a + m) * x;
        num = dnum * (a + b + m);
      } else {
        dnum = m * x;
        num = dnum * (b - m);
      }
      const double d_n = num * inv_den;
      const double dd_n = dnum * inv_den;

      // one step of the modified Lentz method, as for G1 with d2d_n = 0.
      // C_new = 1 + d_n / C, with C C_new = C + d_n
      double C_new = 1.0 + d_n / C;
      double C_C_new = C + d_n;
      if (std::fabs(C_new) < small) {
        C_new = small;
        C_C_new = C * small;
      }
      const double inv_C_C_new = 1.0 / C_C_new;
      const double dlog_C_C_new = (C * dlog_C + dd_n) * inv_C_C_new;
      d2log_C = -(dd_n * dlog_C + d_n * d2log_C) * inv_C_C_new;
      dlog_C = (dd_n - d_n * dlog_C) * inv_C_C_new;
      d2log_C -= dlog_C * dlog_C_C_new;

      // D_new = 1 / (1 + d_n D)
      double D_new = 1.0 + d_n * D;
      if (std::fabs(D_new) < small) {
        D_new = small;
      }
      D_new = 1.0 / D_new;
      const double D_D_new = D * D_new;
      const double dlog_D_D_new = dlog_D - D_D_new * (dd_n + d_n * dlog_D);
      d2log_D = -(dd_n * dlog_D + d_n * d2log_D) * D_D_new;
      dlog_D = -(dd_n + d_n * dlog_D) * D_D_new;
      d2log_D += dlog_D * dlog_D_D_new;

      dlog_h += dlog_C + dlog_D;
      d2log_h += d2log_C + d2log_D;
      C = C_new;
      D = D_new;
      const double delta = C_new * D_new;
      h *= delta;
      if (odd) {
        if (std::fabs(delta - 1.0) < tolerance
            && std::fabs(d2log_h - d2log_h_prev)
                   <= tolerance * std::fabs(d2log_h)) {
          break;
        }
        d2log_h_prev = d2log_h;
      }
    }
    const double dlog_I = log1m_x + digamma_ab - digamma_b + dlog_h;
    const double d2log_Q
        = trigamma_ab - boost::math::trigamma(b, boost_policy_t<>());
    return v_b * h * Q * (dlog_I * dlog_I + d2log_Q + d2log_h);
  } else {
    double log_C_a = 0.0;
    double log_C_b = 0.0;
    double Dv_log_C_a = 0.0;
    double Dv_log_C_b = 0.0;
    const double d_1_a = -(x + d_1) / ap1;
    const double d_1_b = -x / ap1;
    double log_D_a = -d_1_a * D;
    double log_D_b = -d_1_b * D;
    // D_v log D = -D D_v d_1
    const double d_1_aa = 2.0 * (x + d_1) / (ap1 * ap1);
    const double d_1_ab = x / (ap1 * ap1);
    const double Dv_d_1 = d_1_a * v_a + d_1_b * v_b;
    double Dv_log_D_a
        = -(d_1_aa * v_a + d_1_ab * v_b) * D - log_D_a * D * Dv_d_1;
    double Dv_log_D_b = -d_1_ab * v_a * D - log_D_b * D * Dv_d_1;
    double log_h_a = log_D_a;
    double log_h_b = log_D_b;
    double Dv_log_h_a = Dv_log_D_a;
    double Dv_log_h_b = Dv_log_D_b;
    double Dv_log_h_a_prev = Dv_log_h_a;
    double Dv_log_h_b_prev = Dv_log_h_b;
    for (int n = 2; n <= 2 * max_iterations + 1; ++n) {
      const int m = n / 2;
      const bool odd = n % 2 == 1;
      const double a_plus_n = a + n;
      const double inv_den = 1.0 / (a_plus_n * (a_plus_n - 1.0));
      const double den_a = 2.0 * a_plus_n - 1.0;
      double num_b;
      double num;
      double num_a = 0.0;
      double num_aa = 0.0;
      double num_ab = 0.0;
      if (odd) {
        num_b = -(a + m) * x;
        num = num_b * (a + b + m);
        num_a = -(2.0 * (a + m) + b) * x;
        num_aa = -2.0 * x;
        num_ab = -x;
      } else {
        num_b = m * x;
        num = num_b * (b - m);
      }
      const double d_n = num * inv_den;
      const double d_n_a = (num_a - d_n * den_a) * inv_den;
      const double d_n_b = num_b * inv_den;
      const double d_n_aa
          = (num_aa - 2.0 * d_n_a * den_a - 2.0 * d_n) * inv_den;
      const double d_n_ab = (num_ab - d_n_b * den_a) * inv_den;
      const double Dv_d_n = d_n_a * v_a + d_n_b * v_b;
      const double Dv_d_n_a = d_n_aa * v_a + d_n_ab * v_b;
      const double Dv_d_n_b = d_n_ab * v_a;

      // one step of the modified Lentz method. In each half, the Dv_
      // derivatives take their old-step part, the partials update in place,
      // then scale D_v log(C C_new) or D_v log(D D_new).
      // C_new = 1 + d_n / C, with C C_new = C + d_n
      double C_new = 1.0 + d_n / C;
      double C_C_new = C + d_n;
      if (std::fabs(C_new) < small) {
        C_new = small;
        C_C_new = C * small;
      }
      const double inv_C_C_new = 1.0 / C_C_new;
      const double Dv_log_C = log_C_a * v_a + log_C_b * v_b;
      const double Dv_log_C_C_new = (C * Dv_log_C + Dv_d_n) * inv_C_C_new;
      Dv_log_C_a
          = (Dv_d_n_a - Dv_d_n * log_C_a - d_n * Dv_log_C_a) * inv_C_C_new;
      log_C_a = (d_n_a - d_n * log_C_a) * inv_C_C_new;
      Dv_log_C_a -= log_C_a * Dv_log_C_C_new;
      Dv_log_C_b
          = (Dv_d_n_b - Dv_d_n * log_C_b - d_n * Dv_log_C_b) * inv_C_C_new;
      log_C_b = (d_n_b - d_n * log_C_b) * inv_C_C_new;
      Dv_log_C_b -= log_C_b * Dv_log_C_C_new;

      // D_new = 1 / (1 + d_n D)
      double D_new = 1.0 + d_n * D;
      if (std::fabs(D_new) < small) {
        D_new = small;
      }
      D_new = 1.0 / D_new;
      const double D_D_new = D * D_new;
      const double Dv_log_D = log_D_a * v_a + log_D_b * v_b;
      const double Dv_log_D_D_new
          = Dv_log_D - D_D_new * (Dv_d_n + d_n * Dv_log_D);
      Dv_log_D_a = -(Dv_d_n_a + Dv_d_n * log_D_a + d_n * Dv_log_D_a) * D_D_new;
      log_D_a = -(d_n_a + d_n * log_D_a) * D_D_new;
      Dv_log_D_a += log_D_a * Dv_log_D_D_new;
      Dv_log_D_b = -(Dv_d_n_b + Dv_d_n * log_D_b + d_n * Dv_log_D_b) * D_D_new;
      log_D_b = -(d_n_b + d_n * log_D_b) * D_D_new;
      Dv_log_D_b += log_D_b * Dv_log_D_D_new;

      log_h_a += log_C_a + log_D_a;
      log_h_b += log_C_b + log_D_b;
      Dv_log_h_a += Dv_log_C_a + Dv_log_D_a;
      Dv_log_h_b += Dv_log_C_b + Dv_log_D_b;
      C = C_new;
      D = D_new;
      const double delta = C_new * D_new;
      h *= delta;
      if (odd) {
        if (std::fabs(delta - 1.0) < tolerance
            && std::fabs(Dv_log_h_a - Dv_log_h_a_prev)
                   <= tolerance * std::fabs(Dv_log_h_a)
            && std::fabs(Dv_log_h_b - Dv_log_h_b_prev)
                   <= tolerance * std::fabs(Dv_log_h_b)) {
          break;
        }
        Dv_log_h_a_prev = Dv_log_h_a;
        Dv_log_h_b_prev = Dv_log_h_b;
      }
    }
    const double I = h * Q;
    const double log_I_a = log_x - 1.0 / a + digamma_ab - digamma_a + log_h_a;
    const double log_I_b = log1m_x + digamma_ab - digamma_b + log_h_b;
    const double Dv_log_I = log_I_a * v_a + log_I_b * v_b;
    // the Hessian of log Q is that of -log B(a, b) plus diag(1 / a^2, 0)
    const double Dv_log_Q_a
        = (trigamma_ab - boost::math::trigamma(a, boost_policy_t<>())
           + 1.0 / (a * a))
              * v_a
          + trigamma_ab * v_b;
    const double Dv_log_Q_b
        = trigamma_ab * v_a
          + (trigamma_ab - boost::math::trigamma(b, boost_policy_t<>())) * v_b;
    return std::pair{I * (log_I_a * Dv_log_I + Dv_log_Q_a + Dv_log_h_a),
                     I * (log_I_b * Dv_log_I + Dv_log_Q_b + Dv_log_h_b)};
  }
}

/**
 * Hessian-vector product H v of the regularized incomplete beta function
 * I_z(a, b), H its Hessian in (a, b): the reverse pass of
 * (d/da I_z(a, b), d/db I_z(a, b)) with output adjoints v.
 *
 * The branch choice follows grad_reg_inc_beta. Under the symmetry relation
 * I_z(a, b) = 1 - I_{1-z}(b, a), H v = -swap(H' swap(v)), H' the Hessian of
 * I_{1-z}(b, a) in its own arguments, so H_aa = -H'_bb and H_bb = -H'_aa.
 *
 * @tparam grads <code>G1</code> for H_aa v_a, <code>G2</code> for H_bb v_b,
 *   <code>BOTH</code> for H v
 * @param a first shape
 * @param b second shape
 * @param z upper limit, in [0, 1]
 * @param v_a a component of v; not used for <code>G2</code>
 * @param v_b b component of v; not used for <code>G1</code>
 * @param digamma_a digamma(a); not used for <code>G2</code>
 * @param digamma_b digamma(b); not used for <code>G1</code>
 * @param digamma_ab digamma(a + b)
 * @return H_aa v_a for <code>G1</code>, H_bb v_b for <code>G2</code>,
 *   ((H v)_a, (H v)_b) for <code>BOTH</code>; 0 at z = 0 and z = 1
 */
template <GRAD_RIB grads>
inline std::conditional_t<grads == GRAD_RIB::BOTH, std::pair<double, double>,
                          double>
inc_beta_hvp(double a, double b, double z, double v_a, double v_b,
             double digamma_a, double digamma_b, double digamma_ab) {
  if (z == 0.0 || z == 1.0) {
    return {};
  }
  // the symmetry relation swaps which shape is requested
  constexpr GRAD_RIB swapped_grads = swap_grad_rib(grads);
  const double mean = a / (a + b);
  const bool use_power_series
      = (z < mean && b * z <= 1.0 && z <= 0.95)
        || (z >= mean && a * (1.0 - z) <= 1.0 && z >= 0.05);
  if (use_power_series) {
    if (z >= mean) {
      const auto hv = inc_beta_hvp_power_series<swapped_grads>(
          b, a, 1.0 - z, v_b, v_a, digamma_b, digamma_a, digamma_ab);
      if constexpr (grads != GRAD_RIB::BOTH) {
        return -hv;
      } else {
        return {-hv.second, -hv.first};
      }
    } else {
      return inc_beta_hvp_power_series<grads>(a, b, z, v_a, v_b, digamma_a,
                                              digamma_b, digamma_ab);
    }
  } else {
    if (z >= (a - 1.0) / (a + b - 2.0)) {
      const auto hv = inc_beta_hvp_continued_fraction<swapped_grads>(
          b, a, 1.0 - z, v_b, v_a, digamma_b, digamma_a, digamma_ab);
      if constexpr (grads != GRAD_RIB::BOTH) {
        return -hv;
      } else {
        return {-hv.second, -hv.first};
      }
    } else {
      return inc_beta_hvp_continued_fraction<grads>(
          a, b, z, v_a, v_b, digamma_a, digamma_b, digamma_ab);
    }
  }
}

}  // namespace internal

/**
 * The requested gradients of the regularized incomplete beta function
 * <code>inc_beta(a, b, z)</code> with respect to <code>a</code> and
 * <code>b</code>, for <code>var</code> arguments.
 *
 * The values are computed in double by the prim overload. The reverse pass
 * does not differentiate through its series or continued fraction: with
 * output adjoints v, the a and b adjoints are H v for H the Hessian of
 * I_z(a, b) in (a, b), computed directly by internal::inc_beta_hvp, and the
 * z adjoint is the closed form
 *
 *   d/dz d/da I_z(a, b) = inc_beta_ddz(a, b, z) (log z - psi(a) + psi(a+b)),
 *   d/dz d/db I_z(a, b) = inc_beta_ddz(a, b, z) (log1p(-z) - psi(b) +
 *                         psi(a+b)).
 *
 * The digamma arguments only make the values cheap and must equal
 * digamma(a), digamma(b) and digamma(a + b); their adjoints are 0 and the
 * a and b adjoints above are total derivatives. The cross derivative needs
 * the digamma value a one-sided call does not take, so the reverse pass of
 * <code>G1</code> computes digamma(b) when b is a var, and that of
 * <code>G2</code> digamma(a) when a is a var. betaAB is not used.
 *
 * Each argument is a <code>var</code> or arithmetic; mixing with
 * <code>fvar</code> goes to the prim overload.
 *
 * @tparam grads_requested <code>G1</code> for the derivative with respect to
 *   a, <code>G2</code> for the derivative with respect to b,
 *   <code>BOTH</code> for both
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @tparam T4 type of digammaA
 * @tparam T5 type of digammaB
 * @tparam T6 type of digammaSum
 * @tparam T7 type of betaAB
 * @param a a
 * @param b b
 * @param z z
 * @param digammaA the value of <code>digamma(a)</code>; not used for
 *   <code>G2</code>
 * @param digammaB the value of <code>digamma(b)</code>; not used for
 *   <code>G1</code>
 * @param digammaSum the value of <code>digamma(a + b)</code>
 * @param betaAB not used
 * @return d/da I_z(a, b) for <code>G1</code>, d/db for <code>G2</code>, the
 *   pair (d/da, d/db) for <code>BOTH</code>
 */
template <GRAD_RIB grads_requested, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7,
          require_t<bool_constant<internal::is_grad_reg_inc_beta_rev_v<
              T1, T2, T3, T4, T5, T6, T7>>>* = nullptr>
inline auto grad_reg_inc_beta(const T1& a, const T2& b, const T3& z,
                              const T4& digammaA, const T5& digammaB,
                              const T6& digammaSum, const T7& /* betaAB */) {
  const double digamma_a = value_of(digammaA);
  const double digamma_b = value_of(digammaB);
  const double digamma_ab = value_of(digammaSum);
  if constexpr (grads_requested == GRAD_RIB::G1) {
    return make_callback_var(
        grad_reg_inc_beta<grads_requested>(value_of(a), value_of(b),
                                           value_of(z), digamma_a, digamma_b,
                                           digamma_ab, 0.0),
        [a, b, z, digamma_a, digamma_ab](auto& vi) {
          const double a_val = value_of(a);
          const double b_val = value_of(b);
          const double z_val = value_of(z);
          if constexpr (is_var_v<T1> && is_var_v<T2>) {
            // the b adjoint is the cross term H_ab v_a, which needs both
            // sides and digamma(b), which the G1 call does not take
            const auto [hv_a, hv_b] = internal::inc_beta_hvp<GRAD_RIB::BOTH>(
                a_val, b_val, z_val, vi.adj(), 0.0, digamma_a, digamma(b_val),
                digamma_ab);
            a.adj() += hv_a;
            b.adj() += hv_b;
          } else if constexpr (is_var_v<T1>) {
            // H_aa v_a does not use digamma(b)
            a.adj() += internal::inc_beta_hvp<GRAD_RIB::G1>(
                a_val, b_val, z_val, vi.adj(), 0.0, digamma_a, 0.0, digamma_ab);
          } else if constexpr (is_var_v<T2>) {
            b.adj() += internal::inc_beta_hvp<GRAD_RIB::BOTH>(
                           a_val, b_val, z_val, vi.adj(), 0.0, digamma_a,
                           digamma(b_val), digamma_ab)
                           .second;
          }
          if constexpr (is_var_v<T3>) {
            if (z_val > 0.0 && z_val < 1.0) {
              const double log_term = std::log(z_val) - digamma_a + digamma_ab;
              z.adj()
                  += vi.adj() * inc_beta_ddz(a_val, b_val, z_val) * log_term;
            }
          }
        });
  } else if constexpr (grads_requested == GRAD_RIB::G2) {
    return make_callback_var(
        grad_reg_inc_beta<grads_requested>(value_of(a), value_of(b),
                                           value_of(z), digamma_a, digamma_b,
                                           digamma_ab, 0.0),
        [a, b, z, digamma_b, digamma_ab](auto& vi) {
          const double a_val = value_of(a);
          const double b_val = value_of(b);
          const double z_val = value_of(z);
          if constexpr (is_var_v<T1> && is_var_v<T2>) {
            // the a adjoint is the cross term H_ab v_b, which needs both
            // sides and digamma(a), which the G2 call does not take
            const auto [hv_a, hv_b] = internal::inc_beta_hvp<GRAD_RIB::BOTH>(
                a_val, b_val, z_val, 0.0, vi.adj(), digamma(a_val), digamma_b,
                digamma_ab);
            a.adj() += hv_a;
            b.adj() += hv_b;
          } else if constexpr (is_var_v<T1>) {
            a.adj() += internal::inc_beta_hvp<GRAD_RIB::BOTH>(
                           a_val, b_val, z_val, 0.0, vi.adj(), digamma(a_val),
                           digamma_b, digamma_ab)
                           .first;
          } else if constexpr (is_var_v<T2>) {
            // H_bb v_b does not use digamma(a)
            b.adj() += internal::inc_beta_hvp<GRAD_RIB::G2>(
                a_val, b_val, z_val, 0.0, vi.adj(), 0.0, digamma_b, digamma_ab);
          }
          if constexpr (is_var_v<T3>) {
            if (z_val > 0.0 && z_val < 1.0) {
              const double log_term = log1p(-z_val) - digamma_b + digamma_ab;
              z.adj()
                  += vi.adj() * inc_beta_ddz(a_val, b_val, z_val) * log_term;
            }
          }
        });
  } else {  // both
    const auto [g1_val, g2_val] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
        value_of(a), value_of(b), value_of(z), digamma_a, digamma_b, digamma_ab,
        0.0);
    // Both outputs feed one reverse pass, which runs after everything that
    // uses them has added into their adjoints.
    var g1 = g1_val;
    var g2 = g2_val;
    reverse_pass_callback([a, b, z, g1, g2, digamma_a, digamma_b,
                           digamma_ab]() {
      const double a_val = value_of(a);
      const double b_val = value_of(b);
      const double z_val = value_of(z);
      const double v_a = g1.adj();
      const double v_b = g2.adj();
      if constexpr (is_var_v<T1> || is_var_v<T2>) {
        const auto [hv_a, hv_b] = internal::inc_beta_hvp<GRAD_RIB::BOTH>(
            a_val, b_val, z_val, v_a, v_b, digamma_a, digamma_b, digamma_ab);
        if constexpr (is_var_v<T1>) {
          a.adj() += hv_a;
        }
        if constexpr (is_var_v<T2>) {
          b.adj() += hv_b;
        }
      }
      if constexpr (is_var_v<T3>) {
        if (z_val > 0.0 && z_val < 1.0) {
          z.adj() += inc_beta_ddz(a_val, b_val, z_val)
                     * (v_a * (std::log(z_val) - digamma_a + digamma_ab)
                        + v_b * (log1p(-z_val) - digamma_b + digamma_ab));
        }
      }
    });
    return std::pair{g1, g2};
  }
}

}  // namespace math
}  // namespace stan
#endif
