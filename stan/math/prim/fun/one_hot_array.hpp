#ifndef STAN_MATH_PRIM_FUN_ONE_HOT_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_ONE_HOT_ARRAY_HPP

#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return an array with 1 in the k-th position and zero elsewhere.
 *
 * @param K size of the array
 * @param k position of the 1 (indexing from 1)
 * @return An array of size K with all elements initialised to zero
 * and a 1 in the k-th position.
 * @throw std::domain_error if K is not positive, or if k is less than 1 or
 * greater than K.
 */
inline std::vector<double> one_hot_array(int K, int k) {
  static const char* function = "one_hot_array";
  check_positive(function, "size", K);
  check_bounded(function, "k", k, 1, K);

  std::vector<double> v(K, 0);
  v[k - 1] = 1;
  return v;
}

}  // namespace math
}  // namespace stan

#endif
