#ifndef STAN_MATH_REV_META_IS_VAR_HPP
#define STAN_MATH_REV_META_IS_VAR_HPP

#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/rev/core/var.hpp>
#include <type_traits>

namespace stan {

namespace internal {
template <typename T>
struct is_var_impl : std::false_type {};

template <typename T>
struct is_var_impl<math::var_value<T>> : std::true_type {};
}  // namespace internal
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier is a var_value.
 */
template <typename T>
struct is_var<T,
              std::enable_if_t<internal::is_var_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
