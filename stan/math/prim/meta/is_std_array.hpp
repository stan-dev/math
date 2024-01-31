#ifndef STAN_MATH_PRIM_META_IS_STD_ARRAY_HPP
#define STAN_MATH_PRIM_META_IS_STD_ARRAY_HPP

#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <array>
#include <type_traits>

namespace stan {
namespace internal {

/** \ingroup type_trait
 * This underlying implementation is used when the type is not an std array.
 */
template <typename T>
struct is_std_array_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implementation has a static member named value when the
 * template type is an std array.
 */
template <std::size_t Size, typename T>
struct is_std_array_impl<std::array<T, Size>> : std::true_type {};

}  // namespace internal

/** \ingroup type_trait
 * Base implementation for checking if type is std array
 */
template <typename T, typename = void>
struct is_std_array : std::false_type {};

/** \ingroup type_trait
 * Checks if the decayed type of T is a standard array.
 */
template <typename T>
struct is_std_array<
    T, std::enable_if_t<internal::is_std_array_impl<std::decay_t<T>>::value>>
    : std::true_type {};

/** \ingroup type_trait
 * Specialization of scalar_type for array to recursively return the inner
 * scalar type.
 *
 * @tparam T type of standard array
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_array<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

/** \ingroup type_trait
 * Template metaprogram class to compute the type of values stored
 * in a standard array.
 *
 * @tparam T type of elements in standard array.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_std_array<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};
}

#endif
