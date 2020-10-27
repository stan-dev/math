#ifndef STAN_MATH_PRIM_FUN_UNITSPACED_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_UNITSPACED_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array of integers in an ordered sequence.
 *
 * This produces an array from low to high (included).
 *
 * @param low smallest integer
 * @param high largest integer
 * @return An array of size (high - low + 1) with elements linearly spaced
 * between low and high.
 * @throw std::domain_error if high is less than low.
 */
inline std::vector<int> unitspaced_array(int low, int high) {
  static const char* function = "unitspaced_array";
  check_greater_or_equal(function, "high", high, low);

  int K = fabs(high - low + 1);

  Eigen::VectorXi v = Eigen::VectorXi::LinSpaced(K, low, high);
  return {&v[0], &v[0] + K};
}

}  // namespace math
}  // namespace stan

#endif
