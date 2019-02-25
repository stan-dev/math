#ifndef STAN_MATH_PRIM_FUN_SORT_INDICES_DESC_HPP
#define STAN_MATH_PRIM_FUN_SORT_INDICES_DESC_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/meta/index_type.hpp>
#include <stanh/prim/fun/sort_indices.hpp>
#include <algorithm>  // std::sort
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
std::vector<int> sort_indices_desc(const C& xs) {
  return sort_indices<false>(xs);
}

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_FUN_SORT_INDICES_DESC_HPP
#define STAN_MATH_PRIM_FUN_SORT_INDICES_DESC_HPP

#include <stanh/prim/fun/Eigen.hpp>
#include <stanh/prim/meta/index_type.hpp>
#include <stanh/prim/fun/sort_indices.hpp>
#include <algorithm>  // std::sort
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
std::vector<int> sort_indices_desc(const C& xs) {
  return sort_indices<false>(xs);
}

}  // namespace math
}  // namespace stan
#endif
