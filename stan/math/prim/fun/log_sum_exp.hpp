#ifndef STAN_MATH_PRIM_ARR_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_ARR_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

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

}  // namespace math
}  // namespace stan

#endif
