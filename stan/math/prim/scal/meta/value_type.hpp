#ifndef STAN_MATH_PRIM_SCAL_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_SCAL_META_VALUE_TYPE_HPP

#include <type_traits>
#include <stan/math/prim/scal/meta/is_complex.hpp>
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

template<typename T>
struct value_type<T, std::enable_if_t<is_complex<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

/** \ingroup type_trait
 * Helper function for accessing underlying type.
 */
template <typename T>
using value_type_t = typename value_type<T>::type;

}  // namespace stan
#endif
