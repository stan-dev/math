#ifndef STAN_MATH_PRIM_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculates the log sum of exponetials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
 *
   \f[
   \mbox{log\_sum\_exp}(x, y) =
   \begin{cases}
     \ln(\exp(x)+\exp(y)) & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial x} =
   \begin{cases}
     \frac{\exp(x)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log\_sum\_exp}(x, y)}{\partial y} =
   \begin{cases}
     \frac{\exp(y)}{\exp(x)+\exp(y)} & \mbox{if } -\infty\leq x, y \leq \infty
 \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a the first variable
 * @param b the second variable
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> log_sum_exp(const T2& a, const T1& b) {
  if (a == NEGATIVE_INFTY) {
    return b;
  }
  if (a == INFTY && b == INFTY) {
    return INFTY;
  }
  if (a > b) {
    return a + log1p_exp(b - a);
  }
  return b + log1p_exp(a - b);
}

/**
 * Return the log of the sum of the exponentiated values of the specified
 * sequence of values.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @param[in] x array of specified values
 * @return The log of the sum of the exponentiated vector values.
 */
inline double log_sum_exp(const std::vector<double>& x) {
  using std::exp;
  using std::log;
  double max = NEGATIVE_INFTY;
  for (double xx : x) {
    if (xx > max) {
      max = xx;
    }
  }

  double sum = 0.0;
  for (size_t ii = 0; ii < x.size(); ii++) {
    if (x[ii] != NEGATIVE_INFTY) {
      sum += exp(x[ii] - max);
    }
  }

  return max + log(sum);
}

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * or a row vector.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param[in] x Matrix of specified values
 * @return The log of the sum of the exponentiated vector values.
 */
template <int R, int C>
double log_sum_exp(const Eigen::Matrix<double, R, C>& x) {
  if (x.size() == 0) {
    return NEGATIVE_INFTY;
  }

  const double max = x.maxCoeff();
  if (!std::isfinite(max)) {
    return max;
  }
  return max + std::log((x.array() - max).exp().sum());
}

}  // namespace math
}  // namespace stan

#endif
