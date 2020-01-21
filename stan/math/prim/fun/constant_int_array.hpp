#ifndef STAN_MATH_PRIM_FUN_CONSTANT_INT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_CONSTANT_INT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an integer array with elements set to the same constant
 *
 * @param K size of the array
 * @param c constant value
 * @return An array integer of size K with all elements initialised to
 * the same constant.
 * @throw std::domain_error if K is negative.
 */
inline std::vector<int> constant_int_array(int K, int c) {
  check_nonnegative("constant_int_array", "size", K);
  return std::vector<int>(K, c);
}

}  // namespace math
}  // namespace stan

#endif
