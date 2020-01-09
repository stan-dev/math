#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Primary template class for metaprogram to compute the type of
 * values stored in a container.
 *
 * tparam T type of container.
 */
template <typename T, typename = void>
struct value_type {
  using type = typename std::decay_t<T>;
};

/** \ingroup type_trait
 * Specialization for pointers returns the underlying value the pointer is
 * pointing to.
 */
template <typename T>
struct value_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  using type = typename value_type<std::remove_pointer<T>>::type;
};

/** \ingroup type_trait
 * Helper function for accessing underlying type.
 */
template <typename T>
using value_type_t = typename value_type<T>::type;

/** \ingroup type_trait
 * Template metaprogram class to compute the type of values stored
 * in a standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

/** \ingroup type_trait
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};

}  // namespace stan

#endif
