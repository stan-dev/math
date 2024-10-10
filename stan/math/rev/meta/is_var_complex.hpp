#ifndef STAN_MATH_REV_META_IS_VAR_COMPLEX_HPP
#define STAN_MATH_REV_META_IS_VAR_COMPLEX_HPP

#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/rev/core/var.hpp>
#include <type_traits>

namespace stan {

namespace internal {
template <typename T>
struct is_var_complex_impl : std::false_type {};

template <typename T>
struct is_var_complex_impl<math::var_value<std::complex<T>>> : std::true_type {};
}  // namespace internal
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier is a var_value.
 */
template <typename T>
struct is_var_complex<T,
              std::enable_if_t<internal::is_var_complex_impl<std::decay_t<T>>::value>>
    : std::true_type {};

}  // namespace stan
#endif
