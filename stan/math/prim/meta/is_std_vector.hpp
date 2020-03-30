#ifndef STAN_MATH_PRIM_META_IS_STD_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_STD_VECTOR_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace internal {

/** \ingroup type_trait
 * This underlying implimentation is used when the type is not an std vector.
 */
template <typename T>
struct is_std_vector_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implimentation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_std_vector_impl<std::vector<Args...>> : std::true_type {};

}  // namespace internal

/** \ingroup type_trait
 * Base implimentation for checking if type is std vector
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

/** \ingroup type_trait
 * Checks if the decayed type of T is a standard vector.
 */
template <typename T>
struct is_std_vector<
    T, std::enable_if_t<internal::is_std_vector_impl<std::decay_t<T>>::value>>
    : std::true_type {};

/** \ingroup type_trait
 * Specialization of scalar_type for vector to recursively return the inner
 * scalar type.
 *
 * @tparam T type of standard vector
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

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

STAN_ADD_REQUIRE_UNARY(std_vector, is_std_vector, require_std);
STAN_ADD_REQUIRE_CONTAINER(std_vector, is_std_vector, require_std);

}  // namespace stan

#endif
