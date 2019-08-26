#ifndef STAN_MATH_PRIM_SCAL_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_INDEX_TYPE_HPP

#include <type_traits>

namespace stan {
namespace math {

/**
 * Primary template class for the metaprogram to compute the index
 * type of a container.
 *
 * Only the specializations have behavior that can be used, and
 * all implement a typedef <code>type</code> for the type of the
 * index given container <code>T</code>.
 *
 * tparam T type of container.
 */
template <typename T, typename = void>
struct index_type : std::integral_constant<int, 1> {};

}  // namespace math
}  // namespace stan

#endif
