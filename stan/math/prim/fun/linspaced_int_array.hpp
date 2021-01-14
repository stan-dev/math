#ifndef STAN_MATH_PRIM_FUN_LINSPACED_INT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_LINSPACED_INT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array of linearly spaced integers.
 *
 * This produces an array from `low` to `high` (inclusive). If `high - low` is
 * greater or equal to `K - 1`, then the integers are evenly spaced. If it is
 * not possible to get from `low` to `high` with a multiple of an integer,
 * `high` is lowered until this is possible.
 *
 * If `K - 1` is greater than `high - low`, then integers are repeated. For
 * instance, `low, low, low + 1, low + 1, ...`. `high` is lowered until `K - 1`
 * is a multiple of `high - low`
 *
 * For `K = 1`, the array will contain the `high` value. For `K = 0` it returns
 * an empty array.
 *
 * @param K size of the array
 * @param low smallest value
 * @param high largest value
 * @return An array of size K with elements linearly spaced between
 * low and high.
 * @throw std::domain_error if K is negative or high is less than low.
 */
inline std::vector<int> linspaced_int_array(int K, int low, int high) {
  static const char* function = "linspaced_int_array";
  check_nonnegative(function, "size", K);
  check_greater_or_equal(function, "high", high, low);
  if (K == 0) {
    return {};
  }

  Eigen::VectorXi v = Eigen::VectorXi::LinSpaced(K, low, high);
  return {v.data(), v.data() + K};
}

}  // namespace math
}  // namespace stan

#endif
