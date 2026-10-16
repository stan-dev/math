#ifndef STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LOG_CDFS_HPP
#define STAN_MATH_PRIM_PROB_BETA_BINOMIAL_LOG_CDFS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma_diff.hpp>
#include <stan/math/prim/fun/exp.hpp>
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
 * The log cdf and the log ccdf of the beta binomial distribution at one
 * point, and their derivatives with respect to (alpha, beta).
 */
template <typename T>
struct beta_binomial_log_cdfs_t {
  T lcdf;
  T lccdf;
  std::array<T, 2> dlcdf;
  std::array<T, 2> dlccdf;
};

// A sum stops when a bound of its remaining terms is below this fraction of
// the partial sum.
constexpr double bb_cdfs_tol = 1e-17;

template <typename T>
inline double bb_cdfs_max_abs(const std::array<T, 2>& g) {
  return std::max(std::fabs(value_of_rec(g[0])), std::fabs(value_of_rec(g[1])));
}

/**
 * Returns log pmf(k) = log choose(N, k) + lbeta(k + alpha, N - k + beta)
 * - lbeta(alpha, beta), without the cancellation of the two lbeta at large
 * shapes.
 */
template <typename T>
inline T bb_cdfs_log_pmf(int k, int N, const T& alpha, const T& beta) {
  return binomial_coefficient_log(N, k)
         + log_beta_ratio(alpha, beta, T(k), T(N - k),
                          log_beta_ratio_denominator(alpha, beta));
}

/**
 * Writes d log pmf(k) / d(alpha, beta) to g, as differences of digamma
 * functions that do not cancel at large arguments.
 */
template <typename T>
inline void bb_cdfs_dlog_pmf(int k, int N, const T& alpha, const T& beta,
                             std::array<T, 2>& g) {
  const T d_total = digamma_diff(T(alpha + beta), T(N));
  g[0] = digamma_diff(alpha, T(k)) - d_total;
  g[1] = digamma_diff(beta, T(N - k)) - d_total;
}

/**
 * Returns whether A + B i, with i from lo to hi, has the sign sign (or is
 * 0) at both ends; then it has that sign on the whole range. pmf(i + 1) /
 * pmf(i) - 1 has the sign of A + B i, A = N alpha - N + 1 - beta,
 * B = 2 - alpha - beta.
 */
inline bool bb_cdfs_same_sign(double a, double b, double lo, double hi,
                              double sign) {
  return sign * (a + b * lo) >= 0.0 && sign * (a + b * hi) >= 0.0;
}

/**
 * Computes log cdf(n) as the sum of the pmf from k = n down to 0 in log
 * scale, and, if Grad, its gradient as sum pmf(k) d log pmf(k) / cdf(n).
 * All terms are positive. The pmf is unimodal for alpha + beta > 2 and
 * U-shaped for alpha + beta < 2. Where it rises on [0, j], the j + 1
 * remaining terms are each at most pmf(j), and the sum stops when they are
 * negligible.
 */
template <bool Grad, typename T>
inline void bb_cdfs_lower_sum(int n, int N, const T& alpha, const T& beta,
                              T& log_cdf, std::array<T, 2>& dlog) {
  const double a_lin = value_of_rec(N * alpha - N + 1.0 - beta);
  const double b_lin = value_of_rec(2.0 - alpha - beta);
  const double smaller = std::min(value_of_rec(alpha), value_of_rec(beta));
  T log_scale = bb_cdfs_log_pmf(n, N, alpha, beta);
  std::array<T, 2> g{0, 0};
  std::array<T, 2> d{0, 0};
  if (Grad) {
    bb_cdfs_dlog_pmf(n, N, alpha, beta, g);
  }
  T t = 1.0;
  T sum = 0.0;
  constexpr double big = 1e250;
  const double log_big = std::log(big);
  for (int k = n;; --k) {
    sum += t;
    if (Grad) {
      d[0] += t * g[0];
      d[1] += t * g[1];
    }
    if (k == 0) {
      break;
    }
    const double j = k - 1.0;
    // pmf(j) / pmf(j + 1)
    const T q = (j + 1.0) * (N - j - 1.0 + beta) / ((N - j) * (j + alpha));
    if (j == 0.0 || bb_cdfs_same_sign(a_lin, b_lin, 0.0, j - 1.0, 1.0)) {
      // |d log pmf| changes by at most 2 (log1p(N / smaller) + 1 / smaller)
      // over the remaining terms
      const double rest = value_of_rec((j + 1.0) * t * q);
      const double g_bound
          = Grad ? bb_cdfs_max_abs(g)
                       + 2.0 * (std::log1p(N / smaller) + 1.0 / smaller) + 1.0
                 : 1.0;
      if (rest * g_bound <= bb_cdfs_tol * value_of_rec(sum)) {
        break;
      }
    }
    if (Grad) {
      g[0] -= 1.0 / (j + alpha);
      g[1] += 1.0 / (N - j - 1.0 + beta);
    }
    t *= q;
    if (value_of_rec(t) > big || value_of_rec(sum) > big) {
      t /= big;
      sum /= big;
      if (Grad) {
        d[0] /= big;
        d[1] /= big;
      }
      log_scale += log_big;
    }
  }
  log_cdf = log_scale + log(sum);
  if (Grad) {
    dlog[0] = d[0] / sum;
    dlog[1] = d[1] / sum;
  }
}

/**
 * Computes log ccdf(n) as the sum of the pmf from k = n + 1 up to N in log
 * scale, and, if Grad, its gradient. Where the pmf falls on [k + 1, N], the
 * N - k remaining terms are each at most pmf(k + 1), and the sum stops when
 * they are negligible.
 */
template <bool Grad, typename T>
inline void bb_cdfs_upper_sum(int n, int N, const T& alpha, const T& beta,
                              T& log_ccdf, std::array<T, 2>& dlog) {
  const double a_lin = value_of_rec(N * alpha - N + 1.0 - beta);
  const double b_lin = value_of_rec(2.0 - alpha - beta);
  const double smaller = std::min(value_of_rec(alpha), value_of_rec(beta));
  T log_scale = bb_cdfs_log_pmf(n + 1, N, alpha, beta);
  std::array<T, 2> g{0, 0};
  std::array<T, 2> d{0, 0};
  if (Grad) {
    bb_cdfs_dlog_pmf(n + 1, N, alpha, beta, g);
  }
  T t = 1.0;
  T sum = 0.0;
  constexpr double big = 1e250;
  const double log_big = std::log(big);
  for (int k = n + 1;; ++k) {
    sum += t;
    if (Grad) {
      d[0] += t * g[0];
      d[1] += t * g[1];
    }
    if (k == N) {
      break;
    }
    // pmf(k + 1) / pmf(k)
    const T rho = (N - k) * (k + alpha) / ((k + 1.0) * (N - k - 1.0 + beta));
    if (k + 1 == N || bb_cdfs_same_sign(a_lin, b_lin, k + 1.0, N - 1.0, -1.0)) {
      const double rest = value_of_rec((N - k) * t * rho);
      const double g_bound
          = Grad ? bb_cdfs_max_abs(g)
                       + 2.0 * (std::log1p(N / smaller) + 1.0 / smaller) + 1.0
                 : 1.0;
      if (rest * g_bound <= bb_cdfs_tol * value_of_rec(sum)) {
        break;
      }
    }
    if (Grad) {
      g[0] += 1.0 / (k + alpha);
      g[1] -= 1.0 / (N - k - 1.0 + beta);
    }
    t *= rho;
    if (value_of_rec(t) > big || value_of_rec(sum) > big) {
      t /= big;
      sum /= big;
      if (Grad) {
        d[0] /= big;
        d[1] /= big;
      }
      log_scale += log_big;
    }
  }
  log_ccdf = log_scale + log(sum);
  if (Grad) {
    dlog[0] = d[0] / sum;
    dlog[1] = d[1] / sum;
  }
}

/**
 * Returns the log cdf and the log ccdf of the beta binomial distribution at
 * 0 <= n < N, and, if Grad, their derivatives with respect to (alpha,
 * beta). The tail with fewer terms (n + 1 below, N - n above) is summed
 * first. If it is at most 1/2, the other tail is its complement, which then
 * loses no digits; else the other tail is summed too, and if that one is at
 * most 1/2, the first tail is replaced by its complement. A sum near 1
 * whose first term is tiny keeps the rounding error of the log of that term
 * (about 1e-13 for a log of -500); the complement of a small tail does not.
 *
 * @tparam Grad whether to compute the derivatives
 * @tparam T type of the parameters
 * @param n outcome, 0 <= n < N
 * @param N number of trials
 * @param alpha first shape, positive
 * @param beta second shape, positive
 * @return log cdf, log ccdf and, if Grad, their derivatives
 */
template <bool Grad, typename T>
inline beta_binomial_log_cdfs_t<T> beta_binomial_log_cdfs(int n, int N,
                                                          const T& alpha,
                                                          const T& beta) {
  beta_binomial_log_cdfs_t<T> out{0, 0, {0, 0}, {0, 0}};
  const double log_half = -LOG_TWO;
  auto complement
      = [](const T& l, const std::array<T, 2>& dl, T& c, std::array<T, 2>& dc) {
          c = log1m_exp(l);
          if (Grad) {
            const T ratio = exp(l - c);
            dc[0] = -dl[0] * ratio;
            dc[1] = -dl[1] * ratio;
          }
        };
  if (n + 1 <= N - n) {
    bb_cdfs_lower_sum<Grad>(n, N, alpha, beta, out.lcdf, out.dlcdf);
    if (value_of_rec(out.lcdf) <= log_half) {
      complement(out.lcdf, out.dlcdf, out.lccdf, out.dlccdf);
    } else {
      bb_cdfs_upper_sum<Grad>(n, N, alpha, beta, out.lccdf, out.dlccdf);
      if (value_of_rec(out.lccdf) <= log_half) {
        complement(out.lccdf, out.dlccdf, out.lcdf, out.dlcdf);
      }
    }
  } else {
    bb_cdfs_upper_sum<Grad>(n, N, alpha, beta, out.lccdf, out.dlccdf);
    if (value_of_rec(out.lccdf) <= log_half) {
      complement(out.lccdf, out.dlccdf, out.lcdf, out.dlcdf);
    } else {
      bb_cdfs_lower_sum<Grad>(n, N, alpha, beta, out.lcdf, out.dlcdf);
      if (value_of_rec(out.lcdf) <= log_half) {
        complement(out.lcdf, out.dlcdf, out.lccdf, out.dlccdf);
      }
    }
  }
  return out;
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
