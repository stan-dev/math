#ifndef STAN_MATH_PRIM_META_VOID_T_HPP
#define STAN_MATH_PRIM_META_VOID_T_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <type_traits>

namespace stan {

// Dummy struct that has a type void
template <typename... Ts>
struct make_void {
  typedef void type;
};

/**
 * Utility metafunction that maps a sequence of any types to the type void
 * This metafunction is used in template metaprogramming to detect ill-formed
 * types or the validity of an expression in an SFINAE context:
 */
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

}  // namespace stan
#endif
