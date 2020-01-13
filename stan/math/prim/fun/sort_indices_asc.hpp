#ifndef STAN_MATH_PRIM_FUN_SORT_INDICES_ASC_HPP
#define STAN_MATH_PRIM_FUN_SORT_INDICES_ASC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sort_indices.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a sorted copy of the argument container in ascending order.
 *
 * @tparam C type of container
 * @param xs Container to sort
 * @return sorted version of container
 */
template <typename C>
std::vector<int> sort_indices_asc(const C& xs) {
  return sort_indices<true>(xs);
}

}  // namespace math
}  // namespace stan

#endif
