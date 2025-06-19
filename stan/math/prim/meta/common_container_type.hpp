#ifndef STAN_MATH_PRIM_META_COMMON_CONTAINER_TYPE_HPP
#define STAN_MATH_PRIM_META_COMMON_CONTAINER_TYPE_HPP

#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/meta/is_detected.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <type_traits>

namespace stan {
namespace internal {
template <typename T1, typename T2, typename = void, typename = void>
struct common_container_type_impl;

template <typename T1, typename T2>
struct common_container_type_impl<T1, T2, require_stan_scalar_t<T1>,
                                  require_stan_scalar_t<T2>> {
  using type = return_type_t<T1, T2>;
};

template <typename T1, typename T2>
struct common_container_type_impl<T1, T2, require_container_t<T1>,
                                  require_container_t<T2>> {
  using return_t = return_type_t<T1, T2>;
  using container_type_1 = math::promote_scalar_t<return_t, plain_type_t<T1>>;
  using container_type_2 = math::promote_scalar_t<return_t, plain_type_t<T2>>;
  using type =
    std::conditional_t<
      std::is_same<container_type_1, container_type_2>::value,
      container_type_1,
      void
    >;
};

template <typename T1, typename T2>
struct common_container_type_impl<T1, T2, require_stan_scalar_t<T1>,
                                  require_container_t<T2>> {
  using type = math::promote_scalar_t<return_type_t<T1, T2>, plain_type_t<T2>>;
};

template <typename T1, typename T2>
struct common_container_type_impl<T1, T2, require_container_t<T1>,
                                  require_stan_scalar_t<T2>> {
  using type = math::promote_scalar_t<return_type_t<T1, T2>, plain_type_t<T1>>;
};
}

template <typename... Ts>
struct common_container_type;

template <typename T>
struct common_container_type<T> {
  using type = typename internal::common_container_type_impl<
      T, double>::type;  // Use double for base case
};

/**
 * Determine the common container type for a variadic list of types.
 * If all types are scalars, then the common scalar type is returned.
 * If all container types the same, but not necessarily the same scalar type,
 * the common container type with the common scalar type is returned.
 *
 * If different container types are present, the result is `void`.
 */
template <typename T1, typename... Ts>
struct common_container_type<T1, Ts...> {
  using type = typename internal::common_container_type_impl<
      T1, typename common_container_type<Ts...>::type>::type;
};

template <typename... Ts>
using common_container_t = typename common_container_type<Ts...>::type;

}  // namespace stan

#endif  // STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
