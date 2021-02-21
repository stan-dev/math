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
struct is_stan_closure<T, void_t<typename T::captured_scalar_t__>>
    : std::true_type {};

template <typename T>
struct scalar_type<T, void_t<typename T::captured_scalar_t__>> {
  using type = typename T::captured_scalar_t__;
};

STAN_ADD_REQUIRE_UNARY(stan_closure, is_stan_closure, general_types);

template <typename T, typename = void>
struct fn_return_type {
  using type = double;
};

template <typename T>
struct fn_return_type<T, void_t<typename T::captured_scalar_t__>> {
  using type = typename T::captured_scalar_t__;
};

/**
 * Convenience type for the return type of the specified template
 * parameters.
 *
 * @tparam F callable type
 * @tparam Ts sequence of types
 * @see return_type
 * @ingroup type_trait
 */
template <typename F, typename... Args>
using fn_return_type_t
    = return_type_t<typename fn_return_type<F>::type, Args...>;

template <typename T, bool Ref, typename = void>
struct capture_type;

template <typename T>
struct capture_type<T, true> {
  using type = const T&;
};

template <typename T>
struct capture_type<T, false,
                    require_not_stan_closure_t<std::remove_reference_t<T>>> {
  using type = std::remove_reference_t<T>;
};

template <typename T>
struct capture_type<T, false,
                    require_stan_closure_t<std::remove_reference_t<T>>> {
  using type = typename std::remove_reference_t<T>::CopyOf__;
};

/**
 * Type for things captured either by const reference or by copy.
 *
 * @tparam T type of object being captured
 * @tparam Ref true if reference, false if copy
 */
template <typename T, bool Ref>
using capture_type_t = typename capture_type<T, Ref>::type;

}  // namespace stan

#endif
