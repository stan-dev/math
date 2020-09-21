#ifndef STAN_MATH_REV_META_IS_ARENA_MATRIX_HPP
#define STAN_MATH_REV_META_IS_ARENA_MATRIX_HPP

#include <stan/math/rev/meta/arena_type.hpp>
#include <type_traits>

namespace stan {
namespace internal {
template <typename T>
struct is_arena_matrix_impl : std::false_type {};
template <typename T>
struct is_arena_matrix_impl<math::arena_matrix<T>> : std::true_type {};
}  // namespace internal
/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is `arena_matrix<T>`
 */
template <typename T>
struct is_arena_matrix<
    T, require_t<internal::is_arena_matrix_impl<std::decay_t<T>>>>
    : std::true_type {};

}  // namespace stan
#endif
