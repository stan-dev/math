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
 * @tparam T type of container.
 */
template <typename T, typename = void>
struct index_type {};

/**
 * Specialization of index_type for pointers.
 *
 * @tparam T type of container.
 */
template <typename T>
struct index_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  using type = int;
};

template <typename T>
using index_type_t = typename index_type<T>::type;

}  // namespace math
}  // namespace stan

#endif
