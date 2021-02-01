#ifndef STAN_MATH_PRIM_FUN_LINSPACED_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_LINSPACED_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of linearly spaced elements.
 *
 * This produces a vector from low to high (inclusive) with elements spaced
 * as (high - low) / (K - 1). For K=1, the vector will contain the high value;
 * for K=0 it returns an empty vector.
 *
 * @param K size of the vector
 * @param low smallest value
 * @param high largest value
 * @return A vector of size K with elements linearly spaced between
 * low and high.
 * @throw std::domain_error if K is negative, if low is nan or infinite,
 * if high is nan or infinite, or if high is less than low.
 */
inline auto linspaced_vector(int K, double low, double high) {
  static const char* function = "linspaced_vector";
  check_nonnegative(function, "size", K);
  check_finite(function, "low", low);
  check_finite(function, "high", high);
  check_greater_or_equal(function, "high", high, low);

  return Eigen::VectorXd::LinSpaced(K, low, high);
}

}  // namespace math
}  // namespace stan

#endif
