#ifndef STAN_MATH_PRIM_SCAL_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_VALUE_TYPE_HPP

#include <type_traits>

namespace stan {

/**
 * Primary template class for metaprogram to compute the type of
 * values stored in a container.
 *
 * tparam T type of container.
 */
template <typename T, typename = void>
struct value_type {
  using type = typename std::decay_t<T>;
};

/**
 * Specialization for pointers returns the underlying value the pointer is
 * pointing to.
 */
template <typename T>
struct value_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  using type = typename value_type<std::remove_pointer<T>>::type;
};

/**
 * Helper function for accessing underlying type.
 */
template <typename T>
using value_type_t = typename value_type<T>::type;

}  // namespace stan
#endif
