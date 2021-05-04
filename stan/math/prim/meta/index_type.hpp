#ifndef STAN_MATH_PRIM_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_META_INDEX_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/** \ingroup type_trait
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

/** \ingroup type_trait
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

/** \ingroup type_trait
 * Template metaprogram class to compute the type of index for a
 * standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct index_type<T, require_std_vector_t<T>> {
  using type = typename std::decay_t<T>::size_type;
};

/** \ingroup type_trait
 * Template metaprogram defining typedef for the type of index for
 * an Eigen matrix, vector, or row vector.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct index_type<T, require_eigen_t<T>> {
  using type = typename std::decay_t<T>::Index;
};

}  // namespace math

}  // namespace stan

#endif
