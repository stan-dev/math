#ifndef STAN_MATH_REV_META_IS_VAR_HPP
#define STAN_MATH_REV_META_IS_VAR_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier is a var.
 */
template <typename T>
struct is_var<T,
              std::enable_if_t<std::is_same<math::var, std::decay_t<T>>::value>>
    : std::true_type {};

namespace internal {

/** \ingroup type_trait
 * This underlying implementation is used when the type is not an std vector.
 */
template <typename T>
struct is_var_value_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implementation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_var_value_impl<math::var_value<Args...>> : std::true_type {};

}  // namespace internal

template <typename T>
struct is_var_value<T, std::enable_if_t<internal::is_var_value_impl<std::decay_t<T>>::value>> {};

}  // namespace stan
#endif
