#ifndef STAN_MATH_PRIM_FUN_ZEROS_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_ZEROS_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array of zeros.
 *
 * @param K size of the array
 * @return an array of size K with all elements initialized to 0.
 * @throw std::domain_error if K is negative.
 */
inline std::vector<double> zeros_array(int K) {
  check_nonnegative("zeros_array", "size", K);
  return std::vector<double>(K, 0);
}

}  // namespace math
}  // namespace stan

#endif
