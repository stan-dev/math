#ifndef STAN_MATH_PRIM_FUN_ONES_INT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_ONES_INT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an integer array of ones.
 *
 * @param K size of the array
 * @return An integer array of size K with all elements initialised to 1.
 * @throw std::domain_error if K is negative.
 */
inline std::vector<int> ones_int_array(int K) {
  check_nonnegative("ones_int_array", "size", K);
  return std::vector<int>(K, 1);
}

}  // namespace math
}  // namespace stan

#endif
