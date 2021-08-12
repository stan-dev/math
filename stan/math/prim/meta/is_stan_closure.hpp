#ifndef STAN_MATH_PRIM_META_IS_STAN_CLOSURE_HPP
#define STAN_MATH_PRIM_META_IS_STAN_CLOSURE_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/void_t.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>

#include <type_traits>

namespace stan {

/**
 * Checks if type is a closure object.
 * @tparam The type to check
 * @ingroup type_trait
 */
template <typename T, typename = void>
struct is_stan_closure : std::false_type {};

template <typename T>
struct is_stan_closure<T, void_t<typename std::decay_t<T>::return_scalar_t_>>
    : std::true_type {};

STAN_ADD_REQUIRE_UNARY(stan_closure, is_stan_closure, general_types);

template <typename T>
struct scalar_type<T, require_stan_closure_t<T>> {
  using type = typename std::decay_t<T>::return_scalar_t_;
};

template <typename T, bool Ref, typename = void>
struct closure_return_type;

template <typename T>
struct closure_return_type<T, true> {
  using type = const std::decay_t<T>&;
};

template <typename T>
struct closure_return_type<T, false, require_not_stan_closure_t<T>> {
  using type = std::remove_reference_t<T>;
};

template <typename T>
struct closure_return_type<T, false, require_stan_closure_t<T>> {
  using type = typename std::remove_reference_t<T>::Base_;
};

/**
 * Type for things captured either by const reference or by copy.
 *
 * @tparam T type of object being captured
 * @tparam Ref true if reference, false if copy
 */
template <typename T, bool Ref>
using closure_return_type_t = typename closure_return_type<T, Ref>::type;

}  // namespace stan

#endif
