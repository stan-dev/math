#ifndef STAN_MATH_PRIM_FUN_LOG_RISING_FACTORIAL_RATIO_HPP
#define STAN_MATH_PRIM_FUN_LOG_RISING_FACTORIAL_RATIO_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma_stirling_diff.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * Return log1p(t) - t for 0 <= t < 1 without cancellation. For t <= 1/4
 * use log1p(t) = 2 atanh(r), r = t / (2 + t), whose series gives
 * log1p(t) - t = r (2 r^2 (1/3 + r^2/5 + ...) - t) with r^2 <= 1/81.
 *
 * @tparam T type of the argument
 * @param t argument in [0, 1)
 * @return log1p(t) - t
 */
template <typename T>
inline T log1pmx(const T& t) {
  if (t > 0.25) {
    return log1p(t) - t;
  }
  const T r = t / (2.0 + t);
  const T r2 = r * r;
  T s = 1.0 / 21.0;
  s = 1.0 / 19.0 + r2 * s;
  s = 1.0 / 17.0 + r2 * s;
  s = 1.0 / 15.0 + r2 * s;
  s = 1.0 / 13.0 + r2 * s;
  s = 1.0 / 11.0 + r2 * s;
  s = 1.0 / 9.0 + r2 * s;
  s = 1.0 / 7.0 + r2 * s;
  s = 1.0 / 5.0 + r2 * s;
  s = 1.0 / 3.0 + r2 * s;
  return r * (2.0 * r2 * s - t);
}

/**
 * Return the log of the ratio of two rising factorials with the same
 * number of factors,
 *
 *   log((x)_k / (x + b)_k) = lgamma(x + k) - lgamma(x)
 *                            - lgamma(x + b + k) + lgamma(x + b),
 *
 * for x > 0, b >= 0, k >= 0, without the cancellation of the four lgamma
 * terms when x is large. The expression is symmetric in b and k; let
 * m = min(b, k), M = max(b, k) and y = x + M.
 *
 * For x >= lgamma_stirling_diff_useful, write every lgamma through its
 * Stirling form. The terms linear in the arguments cancel exactly and
 *
 *   result = m log1p(-M / (y + m)) + T(x, m) - T(y, m)
 *            + lgamma_stirling_diff(x + m) - lgamma_stirling_diff(x)
 *            - lgamma_stirling_diff(y + m) + lgamma_stirling_diff(y),
 *
 * with T(z, m) = (z - 1/2) log1p(m / z). For m < z, T(z, m) is about m;
 * it is split into m and z (log1p(m / z) - m / z) - log1p(m / z) / 2, and
 * the m parts of the two terms cancel exactly.
 *
 * For x < lgamma_stirling_diff_useful, use
 * lbeta(m, y) - lbeta(m, x), from lgamma(z + m) - lgamma(z)
 * = lgamma(m) - lbeta(m, z).
 *
 * @tparam T type of the arguments
 * @param x first argument, positive
 * @param b shift of the second rising factorial, nonnegative
 * @param k number of factors, nonnegative
 * @return log((x)_k / (x + b)_k)
 */
template <typename T>
inline T log_rising_factorial_ratio(const T& x, const T& b, const T& k) {
  const T m = (b < k) ? b : k;
  const T big = (b < k) ? k : b;
  if (m == 0) {
    return T(0.0);
  }
  const T y = x + big;
  if (x < lgamma_stirling_diff_useful) {
    return lbeta(m, y) - lbeta(m, x);
  }
  const T t_x = m / x;
  const T t_y = m / y;
  T t_diff;
  if (t_x < 1.0) {
    // m < x < y: both split, the two m parts cancel
    t_diff
        = x * log1pmx(t_x) - y * log1pmx(t_y) - 0.5 * (log1p(t_x) - log1p(t_y));
  } else if (t_y < 1.0) {
    // x <= m < y: only the second term is split
    t_diff = (x - 0.5) * log1p(t_x) - m - (y * log1pmx(t_y) - 0.5 * log1p(t_y));
  } else {
    t_diff = (x - 0.5) * log1p(t_x) - (y - 0.5) * log1p(t_y);
  }
  // m log((x + m) / (y + m)); log1p form while the ratio is near 1
  const T shift_frac = big / (y + m);
  const T log_ratio
      = (shift_frac < 0.5) ? log1p(-shift_frac) : log((x + m) / (y + m));
  return m * log_ratio + t_diff + lgamma_stirling_diff(x + m)
         - lgamma_stirling_diff(x) - lgamma_stirling_diff(y + m)
         + lgamma_stirling_diff(y);
}

/**
 * Return (x - 1/2) log1p(k / x) - c, where c is k when k < x and 0
 * otherwise, and add c to `k_removed`. For k < x the term is about k and is
 * formed as x (log1p(k / x) - k / x) - log1p(k / x) / 2, which is of the
 * size k^2 / x. The counts k are integers, so their sum is exact.
 *
 * @tparam T type of the arguments
 * @param x shape, at least lgamma_stirling_diff_useful
 * @param k count, a nonnegative integer
 * @param[in, out] k_removed sum of the removed counts
 * @return the term without the removed count
 */
template <typename T>
inline T stirling_log1p_term(const T& x, const T& k, T& k_removed) {
  if (k == 0) {
    return T(0.0);
  }
  const T t = k / x;
  if (t < 1.0) {
    k_removed += k;
    return x * log1pmx(t) - 0.5 * log1p(t);
  }
  return (x - 0.5) * log1p(t);
}

/**
 * Return an estimate of the size of the terms that `lbeta(a, b)` adds up,
 * for a, b >= lgamma_stirling_diff_useful: min(a, b) (1 + log1p(max / min)).
 * Only used to choose between two forms; the accuracy is not important.
 *
 * @param a first argument
 * @param b second argument
 * @return the size estimate
 */
inline double lbeta_terms_size(double a, double b) {
  const double small = std::fmin(a, b);
  const double large = std::fmax(a, b);
  return small * (1.0 + std::log1p(large / small));
}

/**
 * Return the part of log_beta_ratio() that depends on the shapes only, so
 * that a caller can compute it once per pair of shapes: the three
 * lgamma_stirling_diff terms when both shapes are at least
 * lgamma_stirling_diff_useful, and lbeta(alpha, beta) otherwise.
 *
 * @tparam T type of the arguments
 * @param alpha first shape, positive
 * @param beta second shape, positive
 * @return the shape-only part
 */
template <typename T>
inline T log_beta_ratio_denominator(const T& alpha, const T& beta) {
  if (alpha >= lgamma_stirling_diff_useful
      && beta >= lgamma_stirling_diff_useful) {
    return lgamma_stirling_diff(alpha) + lgamma_stirling_diff(beta)
           - lgamma_stirling_diff(alpha + beta);
  }
  return lbeta(alpha, beta);
}

/**
 * Return lbeta(alpha + n, beta + m) - lbeta(alpha, beta), the log of the
 * ratio of rising factorials (alpha)_n (beta)_m / (alpha + beta)_(n + m),
 * for shapes alpha, beta > 0 and integer counts n, m >= 0.
 *
 * When both shapes are large, each lbeta is of the order of the shapes
 * while the difference is of the order of n + m, so the plain difference
 * keeps no correct digits from shapes near 1e15. For
 * shapes of at least lgamma_stirling_diff_useful, write every lgamma as
 * (x - 1/2) log(x) - x + log(2 pi) / 2 + lgamma_stirling_diff(x). With
 * N = n + m and s = alpha + beta, the terms linear in the arguments cancel
 * exactly and the rest is
 *
 *   T(alpha, n) + T(beta, m) - T(s, N)
 *   + n log((alpha + n) / (s + N)) + m log((beta + m) / (s + N))
 *   + the six lgamma_stirling_diff terms,
 *
 * with T(x, k) = (x - 1/2) log1p(k / x), split by stirling_log1p_term().
 * This form has rounding errors of the size of N; the lbeta form has
 * rounding errors of the size of the smaller shape. The function takes the
 * form with the smaller terms.
 *
 * @tparam T type of the arguments
 * @param alpha first shape, positive
 * @param beta second shape, positive
 * @param n first count, a nonnegative integer
 * @param m second count, a nonnegative integer
 * @param denominator log_beta_ratio_denominator(alpha, beta)
 * @return lbeta(alpha + n, beta + m) - lbeta(alpha, beta)
 */
template <typename T>
inline T log_beta_ratio(const T& alpha, const T& beta, const T& n, const T& m,
                        const T& denominator) {
  const bool large_shapes = alpha >= lgamma_stirling_diff_useful
                            && beta >= lgamma_stirling_diff_useful;
  if (large_shapes) {
    const T total_count = n + m;
    const T alpha_plus_beta = alpha + beta;
    const T total = alpha_plus_beta + total_count;
    T removed_plus(0.0);
    T removed_minus(0.0);
    const T term_alpha = stirling_log1p_term(alpha, n, removed_plus);
    const T term_beta = stirling_log1p_term(beta, m, removed_plus);
    const T term_total
        = stirling_log1p_term(alpha_plus_beta, total_count, removed_minus);
    // n log(p) + m log(q) with p + q = 1: take the log1m form for the larger
    // of p and q, so that a ratio near 1 keeps its digits
    const T p = (alpha + n) / total;
    const T q = (beta + m) / total;
    const T log_n = (p < q) ? n * log(p) : n * log1p(-q);
    const T log_m = (p < q) ? m * log1p(-p) : m * log(q);
    const double size_stirling = std::fabs(value_of_rec(term_alpha))
                                 + std::fabs(value_of_rec(term_beta))
                                 + std::fabs(value_of_rec(term_total))
                                 + std::fabs(value_of_rec(log_n))
                                 + std::fabs(value_of_rec(log_m));
    const double size_lbeta
        = lbeta_terms_size(value_of_rec(alpha), value_of_rec(beta))
          + lbeta_terms_size(value_of_rec(alpha + n), value_of_rec(beta + m));
    if (size_stirling <= size_lbeta) {
      return term_alpha + term_beta - term_total
             + (removed_plus - removed_minus) + log_n + log_m
             + lgamma_stirling_diff(alpha + n) + lgamma_stirling_diff(beta + m)
             - lgamma_stirling_diff(total) - denominator;
    }
    // here denominator holds the Stirling remainders, not lbeta
    return lbeta(alpha + n, beta + m) - lbeta(alpha, beta);
  }
  return lbeta(alpha + n, beta + m) - denominator;
}

/**
 * Return the log pmf of the beta negative binomial distribution at the
 * integer k >= 0,
 *
 *   lgamma(r + k) - lgamma(k + 1) - lgamma(r) + lbeta(alpha + r, beta + k)
 *   - lbeta(alpha, beta),
 *
 * as -log(k) - lbeta(k, beta) + log_rising_factorial_ratio(r, alpha + beta,
 * k) + log_rising_factorial_ratio(alpha, beta, r), which does not cancel
 * when the parameters are large. For k = 0 only the
 * last term remains.
 *
 * @tparam T type of the arguments
 * @param k outcome, a nonnegative integer
 * @param r number of successes, positive
 * @param alpha first shape, positive
 * @param beta second shape, positive
 * @return the log pmf
 */
template <typename T>
inline T beta_neg_binomial_log_pmf(const T& k, const T& r, const T& alpha,
                                   const T& beta) {
  T lp = log_rising_factorial_ratio(alpha, beta, r);
  if (k > 0) {
    lp += log_rising_factorial_ratio(r, T(alpha + beta), k) - log(k)
          - lbeta(k, beta);
  }
  return lp;
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
