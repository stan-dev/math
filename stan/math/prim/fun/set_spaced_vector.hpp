#ifndef STAN_MATH_PRIM_FUN_SET_SPACED_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_SET_SPACED_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector of linearly spaced elements
 *
 * @param K length of the vector
 * @param low smallest value
 * @param high largest value
 * @return A vector of length K with elements linearly spaced between
 * low and high.
 */
Eigen::VectorXd set_spaced_vector(int K, double low, double high) {
  static const char* function = "set_spaced_vector";
  check_nonnegative(function, "length", K);
  check_finite(function, "low", low);
  check_finite(function, "high", high);
  check_greater_or_equal(function, "high", high, low);

  return Eigen::VectorXd::LinSpaced(K, low, high);
}

}  // namespace math
}  // namespace stan

#endif
