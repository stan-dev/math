#ifndef STAN_MATH_PRIM_FUN_RANGE_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_RANGE_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a vector with elements from `low` (included) to `high` (excluded)
 * linearly spaced according to `by`.
 *
 * This produces a sequence of numbers low, low + by, min + 2*by, ...
 * up to but not including high. If low = high or by > high - low, the
 * resulting array will contain only 1 element equal to low.
 *
 * @param low smallest value
 * @param high largest value
 * @param by spacing between elements
 * @return A vector with elements between `low` and `high`, linearly spaced
 * according to `by`.
 * @throw std::domain_error if `low` is nan or infinite, if `high` is nan or
 * infinite, or if `high` is less than `low`, or if `by` is not positive.
 */
inline Eigen::VectorXd range_vector(double low, double high, double by) {
  static const char* function = "range_vector";
  check_finite(function, "low", low);
  check_finite(function, "high", high);
  check_greater_or_equal(function, "high", high, low);
  check_positive_finite(function, "by", by);

  if (low == high) {
    return Eigen::VectorXd::Constant(1, low);
  }

  size_t K = std::ceil((high - low) / by);
  return Eigen::VectorXd::LinSpaced(K, low, low + (K - 1) * by);
}

}  // namespace math
}  // namespace stan

#endif
