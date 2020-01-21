#ifndef STAN_MATH_PRIM_FUN_CONSTANT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_CONSTANT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array with elements set to the same constant
 *
 * @param K size of the array
 * @param c constant value
 * @return An array of size K with all elements initialised to
 * the same constant.
 * @throw std::domain_error if K is negative.
 */
inline std::vector<double> constant_array(int K, double c) {
  check_nonnegative("constant_array", "size", K);
  return std::vector<double>(K, c);
}

}  // namespace math
}  // namespace stan

#endif
