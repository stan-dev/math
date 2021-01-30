#ifndef STAN_MATH_PRIM_META_ARENA_TYPE_HPP
#define STAN_MATH_PRIM_META_ARENA_TYPE_HPP

#include <type_traits>

namespace stan {
namespace math {
// forward declaration
template <typename MatrixType>
class arena_matrix;
}  // namespace math

namespace internal {
template <typename T, typename = void, typename = void>
struct arena_type_impl {};
}  // namespace internal

/**
 * Determines a type that can be used in place of `T` that does any dynamic
 * allocations on the AD stack. This way resulting types are trivially
 * destructible and can be used in vari classes.
 */
template <typename T>
using arena_t = typename internal::arena_type_impl<std::decay_t<T>>::type;

}  // namespace stan
#endif
