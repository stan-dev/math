#ifndef STAN_MATH_PRIM_FUN_LINSPACED_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_LINSPACED_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array of linearly spaced elements.
 *
 * This produces an array from low to high (inclusive) with elements spaced
 * as (high - low) / (K - 1). For K=1, the array will contain the high value;
 * for K=0 it returns an empty array.
 *
 * @param K size of the array
 * @param low smallest value
 * @param high largest value
 * @return An array of size K with elements linearly spaced between
 * low and high.
 * @throw std::domain_error if K is negative, if low is nan or infinite,
 * if high is nan or infinite, or if high is less than low.
 */
inline std::vector<int> linspaced_array(int K, double low, double high) {
  static const char* function = "linspaced_array";
  check_nonnegative(function, "size", K);
  check_finite(function, "low", low);
  check_finite(function, "high", high);
  check_greater_or_equal(function, "high", high, low);

  if (K == 0) {
    return {};
  }

  Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(K, low, high);
  return {&v[0], &v[0] + K};
}

}  // namespace math
}  // namespace stan

#endif
