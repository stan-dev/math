#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

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
