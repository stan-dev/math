#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_TAIL_BOUND_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_TAIL_BOUND_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Return the smallest value of |p + j| over the integers j in [lo, hi].
 *
 * @param p parameter
 * @param lo smallest j
 * @param hi largest j, at least lo
 * @return smallest |p + j|
 */
inline double hypergeometric_3F2_min_abs(double p, double lo, double hi) {
  if ((p + lo) * (p + hi) > 0.0) {
    return std::fmin(std::fabs(p + lo), std::fabs(p + hi));
  }
  // p + j has a zero in [lo, hi]: the smallest value is at one of the two
  // integers next to -p
  const double j = std::floor(-p);
  double min_abs = std::numeric_limits<double>::infinity();
  for (double i : {j, j + 1}) {
    if (i >= lo && i <= hi) {
      min_abs = std::fmin(min_abs, std::fabs(p + i));
    }
  }
  return min_abs;
}

/**
 * Return an upper bound of |x + j| / |y + j| over the integers j in
 * [lo, hi].
 *
 * If neither x + j nor y + j changes sign in [lo, hi], the function of j is
 * a Moebius transformation without a pole there, so it is monotone and has
 * its largest value at lo or hi. Otherwise the bound is the largest
 * |x + j| (at lo or hi) divided by the smallest |y + j|.
 *
 * @param x numerator parameter
 * @param y denominator parameter
 * @param lo smallest j
 * @param hi largest j, at least lo
 * @return upper bound
 */
inline double hypergeometric_3F2_factor_bound(double x, double y, double lo,
                                              double hi) {
  if ((x + lo) * (x + hi) >= 0.0 && (y + lo) * (y + hi) > 0.0) {
    return std::fmax(std::fabs((x + lo) / (y + lo)),
                     std::fabs((x + hi) / (y + hi)));
  }
  return std::fmax(std::fabs(x + lo), std::fabs(x + hi))
         / hypergeometric_3F2_min_abs(y, lo, hi);
}

/**
 * Return an upper bound of the absolute value of the term ratio
 * r_j = z (a1 + j)(a2 + j)(a3 + j) / ((b1 + j)(b2 + j)(1 + j))
 * of the 3F2 series over the integers j in [lo, hi].
 *
 * The ratio is the product of three factors |a_i + j| / |d_l + j|, with
 * d = (b1, b2, 1). Every assignment of the numerator parameters to the
 * denominator parameters gives a bound; the function returns the smallest
 * of the six.
 *
 * @param a numerator parameters
 * @param b denominator parameters
 * @param abs_z absolute value of the argument z
 * @param lo smallest j
 * @param hi largest j, at least lo
 * @return upper bound of |r_j|
 */
inline double hypergeometric_3F2_ratio_bound(const std::array<double, 3>& a,
                                             const std::array<double, 2>& b,
                                             double abs_z, double lo,
                                             double hi) {
  const std::array<double, 3> d{b[0], b[1], 1.0};
  std::array<std::array<double, 3>, 3> factor;
  for (int i = 0; i < 3; ++i) {
    for (int l = 0; l < 3; ++l) {
      factor[i][l] = hypergeometric_3F2_factor_bound(a[i], d[l], lo, hi);
    }
  }
  static constexpr int pairings[6][3]
      = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  double bound = std::numeric_limits<double>::infinity();
  for (const auto& p : pairings) {
    bound
        = std::fmin(bound, factor[0][p[0]] * factor[1][p[1]] * factor[2][p[2]]);
  }
  return abs_z * bound;
}

/**
 * Return upper bounds of sum_{d = 1}^{n} rho^d and sum_{d = 1}^{n} d rho^d.
 *
 * With these sums, a term t of the series and a bound rho of the next n
 * term ratios bound the sum of the next n terms by |t| times the first
 * sum. A derivative term t s_j with |s_j| <= |s| + (j - k) c has the bound
 * |t| (|s| times the first sum + c times the second sum).
 *
 * @param rho bound of the absolute term ratios
 * @param n number of remaining terms
 * @return the two bounds
 */
inline std::pair<double, double> hypergeometric_3F2_tail_sums(double rho,
                                                              double n) {
  if (n <= 0) {
    return {0.0, 0.0};
  }
  if (rho < 1.0) {
    return {
        std::fmin(n * rho, rho / (1.0 - rho)),
        std::fmin(0.5 * n * (n + 1) * rho, rho / ((1.0 - rho) * (1.0 - rho)))};
  }
  const double rho_n = std::pow(rho, n);
  return {n * rho_n, n * n * rho_n};
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
