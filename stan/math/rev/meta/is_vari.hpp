#ifndef STAN_MATH_REV_META_IS_VARI_HPP
#define STAN_MATH_REV_META_IS_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {

namespace internal {
template <typename T>
struct is_vari_impl : std::false_type {};

template <typename T>
struct is_vari_impl<math::vari_value<T>> : std::true_type {};

template <typename T>
struct is_vari_impl<math::vari_view<T>> : std::true_type {};

}  // namespace internal
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a
 * vari_value.
 */
template <typename T>
struct is_vari<T, require_t<internal::is_vari_impl<std::decay_t<T>>>>
    : std::true_type {};

}  // namespace stan
#endif
