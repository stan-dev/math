#ifndef STAN_MATH_PRIM_META_ARENA_TYPE_HPP
#define STAN_MATH_PRIM_META_ARENA_TYPE_HPP

namespace stan {
namespace math {
namespace internal {
template <typename T, typename = void, typename = void>
struct arena_type_impl;
}

/**
 * Determines a type that can be used in place of `T` that does any dynamic
 * allocations on the AD stack. This way resulting types are trivially
 * destructible and can be used in vari classes.
 */
template <typename T>
using arena_t = typename internal::arena_type_impl<std::decay_t<T>, void, void>::type;


template <typename T>
using require_arena_t = require_t<std::is_same<arena_t<std::decay_t<T>>, std::decay_t<T>>>;

template <typename T>
using require_not_arena_t = require_not_t<std::is_same<arena_t<std::decay_t<T>>, std::decay_t<T>>>;

}
}

#endif
