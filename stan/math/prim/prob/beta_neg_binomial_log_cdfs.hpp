#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LOG_CDFS_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_LOG_CDFS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>
#include <stan/math/prim/fun/log_rising_factorial_ratio.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <algorithm>
#include <array>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * The log cdf and the log ccdf of the beta negative binomial distribution at
 * one point, and their derivatives with respect to (r, alpha, beta).
 */
template <typename T>
struct beta_neg_binomial_log_cdfs_t {
  T lcdf;
  T lccdf;
  std::array<T, 3> dlcdf;
  std::array<T, 3> dlccdf;
};

// A sum stops when a bound of its remaining terms is below this fraction of
// the partial sum.
constexpr double bnb_cdfs_tol = 1e-17;
// The complement 1 - p is used for a probability p of at least this size;
// it multiplies the relative error of p by at most (1 - p) / p < 100.
constexpr double bnb_cdfs_complement_min = 1e-2;
constexpr int bnb_cdfs_short_tail_steps = 2000;
// Up to this n the lower sum, with at most n + 1 terms, comes first.
constexpr int bnb_cdfs_lower_first_max_n = 100;
constexpr int bnb_cdfs_thomae_steps = 5000;
constexpr double bnb_cdfs_thomae_max_cancellation = 100;

template <typename T>
inline double bnb_cdfs_max_abs(const std::array<T, 3>& g) {
  return std::max({std::fabs(value_of_rec(g[0])), std::fabs(value_of_rec(g[1])),
                   std::fabs(value_of_rec(g[2]))});
}

/**
 * Writes d log pmf(k) / d(r, alpha, beta) to g, as differences of digamma
 * functions that do not cancel at large arguments.
 */
template <typename T>
inline void bnb_cdfs_dlog_pmf(const T& k, const T& r, const T& alpha,
                              const T& beta, std::array<T, 3>& g) {
  const T dab = digamma_diff(T(alpha + beta), T(k + r));
  g[0] = digamma_diff(r, k) - digamma_diff(T(r + alpha), T(k + beta));
  g[1] = digamma_diff(alpha, r) - dab;
  g[2] = digamma_diff(beta, k) - dab;
}

/**
 * Computes log cdf(n) as the sum of the pmf from k = n down to k = 0 in log
 * scale, and its gradient as sum pmf(k) d log pmf(k) / cdf(n). All terms are
 * positive. pmf(k + 1) / pmf(k) - 1 has the sign of
 * r beta - s - k (1 + alpha), s = r + alpha + beta, so the pmf is unimodal:
 * below the mode the terms fall as k falls, the j remaining terms are each
 * below the current one, and the sum stops when they are negligible.
 */
template <bool Grad, typename T>
inline void bnb_cdfs_lower_sum(int n, const T& r, const T& alpha, const T& beta,
                               int max_steps, const char* function, T& log_cdf,
                               std::array<T, 3>& dlog) {
  const T s = r + alpha + beta;
  const double k_peak = value_of_rec((r * beta - s) / (1.0 + alpha));
  const double smaller = std::min(value_of_rec(r), value_of_rec(beta));
  T log_scale = beta_neg_binomial_log_pmf(T(n), r, alpha, beta);
  std::array<T, 3> g{0, 0, 0};
  std::array<T, 3> d{0, 0, 0};
  if (Grad) {
    bnb_cdfs_dlog_pmf(T(n), r, alpha, beta, g);
  }
  T t = 1.0;
  T sum = 0.0;
  constexpr double big = 1e250;
  const double log_big = std::log(big);
  for (int k = n, steps = 0;; --k, ++steps) {
    if (steps > max_steps) {
      throw_domain_error(function, "number of terms of the lower sum",
                         max_steps, "exceeded ", "");
    }
    sum += t;
    if (Grad) {
      for (int i = 0; i < 3; ++i) {
        d[i] += t * g[i];
      }
    }
    if (k == 0) {
      break;
    }
    const double j = k - 1.0;
    // pmf(j) / pmf(j + 1)
    const T q = (j + 1.0) * (j + s) / ((j + r) * (j + beta));
    if (j < k_peak) {
      // the j + 1 terms below are each below t q; on the way down
      // |d log pmf| changes by at most 3 (log1p(j / smaller) + 1 / smaller)
      const double rest = value_of_rec((j + 1.0) * t * q);
      const double g_bound
          = Grad ? bnb_cdfs_max_abs(g)
                       + 3.0 * (std::log1p(j / smaller) + 1.0 / smaller) + 1.0
                 : 1.0;
      if (rest * g_bound <= bnb_cdfs_tol * value_of_rec(sum)) {
        break;
      }
    }
    if (Grad) {
      const T js = j + s;
      g[0] -= (alpha + beta) / ((j + r) * js);
      g[1] += 1.0 / js;
      g[2] -= (r + alpha) / ((j + beta) * js);
    }
    t *= q;
    if (value_of_rec(t) > big || value_of_rec(sum) > big) {
      t /= big;
      sum /= big;
      if (Grad) {
        for (int i = 0; i < 3; ++i) {
          d[i] /= big;
        }
      }
      log_scale += log_big;
    }
  }
  log_cdf = log_scale + log(sum);
  if (Grad) {
    for (int i = 0; i < 3; ++i) {
      dlog[i] = d[i] / sum;
    }
  }
}

/**
 * Computes log ccdf(n) as the sum of the pmf over k > n, and its gradient.
 * Past the peak of the terms the rest is about t / c for a geometric fall,
 * c = 1 - pmf(k + 1) / pmf(k), and about t (k + s) / alpha for a fall like
 * k^(-1 - alpha); the larger estimate is used. Returns false if the sum
 * does not stop within max_terms terms.
 */
template <bool Grad, typename T>
inline bool bnb_cdfs_tail_sum(int n, const T& r, const T& alpha, const T& beta,
                              int max_terms, T& log_ccdf,
                              std::array<T, 3>& dlog) {
  const T s = r + alpha + beta;
  const double k_peak = value_of_rec((r * beta - s) / (1.0 + alpha));
  const double alpha_dbl = value_of_rec(alpha);
  double k = n + 1.0;
  const T log_p = beta_neg_binomial_log_pmf(T(k), r, alpha, beta);
  std::array<T, 3> g{0, 0, 0};
  std::array<T, 3> d{0, 0, 0};
  if (Grad) {
    bnb_cdfs_dlog_pmf(T(k), r, alpha, beta, g);
  }
  T t = 1.0;
  T sum = 0.0;
  for (int terms = 1; terms <= max_terms; ++terms, k += 1.0) {
    sum += t;
    if (Grad) {
      for (int i = 0; i < 3; ++i) {
        d[i] += t * g[i];
      }
    }
    const T ks = k + s;
    if (k > k_peak) {
      const double c
          = value_of_rec((k * (1.0 + alpha) + s - r * beta) / ((k + 1.0) * ks));
      const double rest
          = value_of_rec(t) * std::max(1.0 / c, value_of_rec(ks) / alpha_dbl);
      const double g_bound
          = Grad ? bnb_cdfs_max_abs(g) + 1.0 + 1.0 / alpha_dbl : 1.0;
      if (rest * g_bound <= bnb_cdfs_tol * value_of_rec(sum)) {
        log_ccdf = log_p + log(sum);
        if (Grad) {
          for (int i = 0; i < 3; ++i) {
            dlog[i] = d[i] / sum;
          }
        }
        return true;
      }
    }
    t *= (k + r) * (k + beta) / ((k + 1.0) * ks);
    if (Grad) {
      g[0] += (alpha + beta) / ((k + r) * ks);
      g[1] -= 1.0 / ks;
      g[2] += (r + alpha) / ((k + beta) * ks);
    }
  }
  return false;
}

/**
 * Returns whether bnb_cdfs_thomae() takes x = r (else x = beta): the smaller
 * of x (alpha + y) gives the smaller first terms.
 */
inline bool bnb_cdfs_thomae_x_is_r(double r, double alpha, double beta) {
  return r * (alpha + beta) <= beta * (alpha + r);
}

/**
 * Returns an estimate of the number of terms of bnb_cdfs_tail_sum() for a
 * geometric fall: log(1 / bnb_cdfs_tol) / c0, c0 = 1 - pmf(n + 2) /
 * pmf(n + 1), plus the steps up to the peak of the terms. The estimate is
 * optimistic where the terms later fall like k^(-1 - alpha); it only
 * selects the method.
 */
inline double bnb_cdfs_tail_sum_terms(int n, double r, double alpha,
                                      double beta) {
  const double s = r + alpha + beta;
  const double k_peak = (r * beta - s) / (1.0 + alpha);
  const double k = n + 1.0;
  const double log_tol = -std::log(bnb_cdfs_tol);
  if (k <= k_peak) {
    return k_peak - k + log_tol;
  }
  const double c0 = (k * (1.0 + alpha) + s - r * beta) / ((k + 1.0) * (k + s));
  return (c0 > 0.0) ? log_tol / c0 : INFTY;
}

/**
 * Returns an estimate of the number of terms of bnb_cdfs_thomae(), or
 * infinity if it exceeds bnb_cdfs_thomae_steps. Over K terms the terms of
 * the series fall by about x log(K) + log1p(K / alpha) + (n + 1) log1p(K /
 * A); the series stops when the fall reaches log(1 / bnb_cdfs_tol) plus the
 * log of the factor max(1, A / (x + n + 1)) of the rest estimate. The
 * estimate only selects the method.
 */
inline double bnb_cdfs_thomae_terms(int n, double r, double alpha,
                                    double beta) {
  const bool x_is_r = bnb_cdfs_thomae_x_is_r(r, alpha, beta);
  const double x = x_is_r ? r : beta;
  const double y = x_is_r ? beta : r;
  const double n1 = n + 1.0;
  const double a = alpha + y + n1;
  const double needed
      = -std::log(bnb_cdfs_tol) + std::log(std::max(1.0, a / (x + n1)));
  for (double k = 16; k <= bnb_cdfs_thomae_steps; k *= 4) {
    const double fall
        = x * std::log(k) + std::log1p(k / alpha) + n1 * std::log1p(k / a);
    if (fall >= needed) {
      return k;
    }
  }
  return INFTY;
}

/**
 * Computes log ccdf(n) and its gradient by a Thomae transformation of the
 * series in ccdf = pmf(n + 1) F, F = 3F2(1, beta + n + 1, r + n + 1;
 * n + 2, alpha + r + beta + n + 1; 1), whose terms fall only like
 * k^(-1 - alpha). With x one of r and beta, y the other one, and
 * A = alpha + y + n + 1:
 *
 *   F = (n + 1) / alpha * (A)_x / (n + 1)_x
 *       * 3F2(1 - x, alpha + y, alpha; alpha + 1, A; 1),
 *
 * where the new series has excess x + n + 1. For x > 1 its terms change
 * sign, so the function returns false if the sum of their absolute values
 * exceeds the sum by more than bnb_cdfs_thomae_max_cancellation, or if the
 * series does not stop within bnb_cdfs_thomae_steps terms. For an integer
 * x the factor 1 - x + k vanishes at k = x - 1: the later terms are 0, but
 * their derivatives with respect to x are not.
 */
template <bool Grad, typename T>
inline bool bnb_cdfs_thomae(int n, const T& r, const T& alpha, const T& beta,
                            T& log_ccdf, std::array<T, 3>& dlog) {
  const bool x_is_r = bnb_cdfs_thomae_x_is_r(
      value_of_rec(r), value_of_rec(alpha), value_of_rec(beta));
  const T& x = x_is_r ? r : beta;
  const T& y = x_is_r ? beta : r;
  const double n1 = n + 1.0;
  const T a = alpha + y + n1;
  const double x_plus_n1 = value_of_rec(x) + n1;
  T t = 1.0;
  T sum = 0.0;
  double abs_sum = 0.0;
  // derivatives of the sum and of log t with respect to (x, y, alpha)
  std::array<T, 3> dsum{0, 0, 0};
  std::array<T, 3> h{0, 0, 0};
  bool polynomial = false;
  T u = 0.0;
  bool done = false;
  for (int k = 0; k < bnb_cdfs_thomae_steps; ++k) {
    const T a1 = 1.0 - x + k;
    const T a2 = alpha + y + k;
    const T a3 = alpha + k;
    const T b1 = alpha + 1.0 + k;
    const T b2 = a + k;
    const T others = a2 * a3 / (b1 * b2 * (k + 1.0));
    if (polynomial) {
      // d t_k / dx = -u_k after the zero factor
      dsum[0] -= u;
      const double abs_q = std::fabs(value_of_rec(a1 * others));
      if (abs_q < 1.0) {
        const double rest
            = std::fabs(value_of_rec(u)) * abs_q
              * std::max(1.0 / (1.0 - abs_q), value_of_rec(b2) / x_plus_n1);
        if (rest <= bnb_cdfs_tol
                        * (std::fabs(value_of_rec(sum))
                           + std::fabs(value_of_rec(dsum[0])))) {
          done = true;
          break;
        }
      }
      u *= a1 * others;
      continue;
    }
    sum += t;
    abs_sum += std::fabs(value_of_rec(t));
    if (Grad) {
      for (int i = 0; i < 3; ++i) {
        dsum[i] += t * h[i];
      }
    }
    if (value_of_rec(a1) == 0.0) {
      if (!Grad) {
        done = true;
        break;
      }
      polynomial = true;
      u = t * others;
      continue;
    }
    const T q = a1 * others;
    const double abs_q = std::fabs(value_of_rec(q));
    if (value_of_rec(a1) > 0.0 && abs_q < 1.0) {
      const double rest
          = std::fabs(value_of_rec(t)) * abs_q
            * std::max(1.0 / (1.0 - abs_q), value_of_rec(b2) / x_plus_n1);
      const double abs_partial = std::fabs(value_of_rec(sum));
      // no sign changes remain: the sum cannot grow by more than rest
      if (abs_sum > bnb_cdfs_thomae_max_cancellation * (abs_partial + rest)) {
        return false;
      }
      const double h_bound = Grad ? bnb_cdfs_max_abs(h) + 1.0 : 1.0;
      if (rest * h_bound <= bnb_cdfs_tol * abs_partial) {
        done = true;
        break;
      }
    }
    if (Grad) {
      h[0] -= 1.0 / a1;
      h[1] += n1 / (a2 * b2);
      h[2] += n1 / (a2 * b2) + 1.0 / (a3 * b1);
    }
    t *= q;
  }
  if (!done || !(value_of_rec(sum) > 0.0)
      || abs_sum > bnb_cdfs_thomae_max_cancellation * value_of_rec(sum)) {
    return false;
  }
  log_ccdf = beta_neg_binomial_log_pmf(T(n1), r, alpha, beta) + std::log(n1)
             - log_rising_factorial_ratio(T(n1), T(alpha + y), x) - log(alpha)
             + log(sum);
  if (Grad) {
    std::array<T, 3> g;
    bnb_cdfs_dlog_pmf(T(n1), r, alpha, beta, g);
    const T d_y_prefactor = digamma_diff(a, x);
    const T d_x = digamma_diff(T(x + n1), T(alpha + y)) + dsum[0] / sum;
    const T d_y = d_y_prefactor + dsum[1] / sum;
    dlog[0] = g[0] + (x_is_r ? d_x : d_y);
    dlog[1] = g[1] + d_y_prefactor - 1.0 / alpha + dsum[2] / sum;
    dlog[2] = g[2] + (x_is_r ? d_y : d_x);
  }
  return true;
}

/**
 * Returns the log cdf and the log ccdf of the beta negative binomial
 * distribution at n >= 0, and, if Grad, their derivatives with respect to
 * (r, alpha, beta). One of the two tails is summed directly and the other
 * one is its complement, so that no result loses digits to cancellation:
 *
 *   - past the mode and for n > bnb_cdfs_lower_first_max_n, the upper tail
 *     first: a short sum of the pmf over k > n, or a Thomae transformation
 *     of its series (fast for large n), the one with the smaller estimate
 *     of its number of terms first;
 *   - else the lower tail, the sum of the pmf from n down to 0;
 *   - the complement where the computed tail is at least 1/2, or the other
 *     tail is at least bnb_cdfs_complement_min;
 *   - else the upper tail by the short sum, the Thomae series or a long sum.
 *
 * @tparam Grad whether to compute the derivatives
 * @tparam T type of the parameters
 * @param n outcome, at least 0
 * @param r number of successes, positive
 * @param alpha prior success, positive
 * @param beta prior failure, positive
 * @param max_steps largest number of terms of a sum
 * @param function name of the calling function, for the error message
 * @return log cdf, log ccdf and, if Grad, their derivatives
 * @throw std::domain_error if a sum needs more than max_steps terms
 */
template <bool Grad, typename T>
inline beta_neg_binomial_log_cdfs_t<T> beta_neg_binomial_log_cdfs(
    int n, const T& r, const T& alpha, const T& beta, int max_steps,
    const char* function) {
  beta_neg_binomial_log_cdfs_t<T> out{0, 0, {0, 0, 0}, {0, 0, 0}};
  const double k_peak
      = value_of_rec((r * beta - (r + alpha + beta)) / (1.0 + alpha));
  const double log_half = -LOG_TWO;
  const double log_complement_min = std::log(bnb_cdfs_complement_min);
  T lc = 0.0;
  std::array<T, 3> d{0, 0, 0};
  bool have_upper = false;
  auto from_upper = [&]() {
    out.lccdf = lc;
    out.lcdf = log1m_exp(lc);
    if (Grad) {
      const T ratio = exp(lc - out.lcdf);
      for (int i = 0; i < 3; ++i) {
        out.dlccdf[i] = d[i];
        out.dlcdf[i] = -d[i] * ratio;
      }
    }
    return out;
  };
  // the fast upper-tail methods, the one with the smaller estimate of its
  // number of terms first; the tail-sum estimate is optimistic for a fall
  // like k^(-1 - alpha), hence the factor 4. A method whose estimate exceeds
  // its budget is not tried.
  auto upper_fast = [&]() {
    const double r_dbl = value_of_rec(r);
    const double alpha_dbl = value_of_rec(alpha);
    const double beta_dbl = value_of_rec(beta);
    const double tail_terms
        = bnb_cdfs_tail_sum_terms(n, r_dbl, alpha_dbl, beta_dbl);
    const double thomae_terms
        = bnb_cdfs_thomae_terms(n, r_dbl, alpha_dbl, beta_dbl);
    auto tail = [&]() {
      return tail_terms <= bnb_cdfs_short_tail_steps
             && bnb_cdfs_tail_sum<Grad>(n, r, alpha, beta,
                                        bnb_cdfs_short_tail_steps, lc, d);
    };
    auto thomae = [&]() {
      return thomae_terms <= bnb_cdfs_thomae_steps
             && bnb_cdfs_thomae<Grad>(n, r, alpha, beta, lc, d);
    };
    if (thomae_terms < 4.0 * tail_terms) {
      return thomae() || tail();
    }
    return tail() || thomae();
  };
  // the lower sum has at most n + 1 terms, and fewer below the mode
  const bool tried_fast = n > k_peak && n > bnb_cdfs_lower_first_max_n;
  if (tried_fast) {
    have_upper = upper_fast();
    // the complement of the upper tail is accurate if the lower tail is
    // not small
    if (have_upper
        && (value_of_rec(lc) <= log_half
            || value_of_rec(log1m_exp(lc)) >= log_complement_min)) {
      return from_upper();
    }
  }
  bnb_cdfs_lower_sum<Grad>(n, r, alpha, beta, max_steps, function, out.lcdf,
                           out.dlcdf);
  const double lcdf_dbl = value_of_rec(out.lcdf);
  if (lcdf_dbl < log_half || -std::expm1(lcdf_dbl) >= bnb_cdfs_complement_min) {
    out.lccdf = log1m_exp(out.lcdf);
    if (Grad) {
      const T ratio = exp(out.lcdf - out.lccdf);
      for (int i = 0; i < 3; ++i) {
        out.dlccdf[i] = -out.dlcdf[i] * ratio;
      }
    }
    return out;
  }
  // the upper tail is below bnb_cdfs_complement_min
  if (!have_upper && !tried_fast) {
    have_upper = upper_fast();
  }
  if (!have_upper) {
    have_upper = bnb_cdfs_tail_sum<Grad>(n, r, alpha, beta, max_steps, lc, d);
  }
  if (!have_upper) {
    throw_domain_error(function, "number of terms of the upper sum", max_steps,
                       "exceeded ", "");
  }
  return from_upper();
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
