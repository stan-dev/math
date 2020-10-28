#ifndef STAN_MATH_PRIM_META_IS_ARENA_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_ARENA_MATRIX_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Defines a static member named value which is defined to be true
 * if the type is `arena_matrix<T>`
 */
template <typename T, typename = void>
struct is_arena_matrix : std::false_type {};

STAN_ADD_REQUIRE_UNARY(arena_matrix, is_arena_matrix, require_eigens_types);
STAN_ADD_REQUIRE_CONTAINER(arena_matrix, is_arena_matrix, require_eigens_types);
STAN_ADD_REQUIRE_UNARY_INNER(arena_matrix, is_arena_matrix,
                             require_eigens_types);

}  // namespace stan
#endif
