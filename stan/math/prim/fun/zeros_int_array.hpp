#ifndef STAN_MATH_PRIM_FUN_ZEROS_INT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_ZEROS_INT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an integer array of zeros.
 *
 * @param K size of the array
 * @return an integer array of size K with all elements initialized to 0.
 * @throw std::domain_error if K is negative.
 */
inline std::vector<int> zeros_int_array(int K) {
  check_nonnegative("zeros_int_array", "size", K);
  return std::vector<int>(K, 0);
}

}  // namespace math
}  // namespace stan

#endif
