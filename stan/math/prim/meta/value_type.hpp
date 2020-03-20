#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Primary template class for metaprogram to compute the type of
 * values stored in a container.
 *
 * tparam T type of container.
 */
template <typename T, typename = void>
struct value_type {
  using type = typename std::decay_t<T>;
};

/** \ingroup type_trait
 * Specialization for pointers returns the underlying value the pointer is
 * pointing to.
 */
template <typename T>
struct value_type<T, std::enable_if_t<std::is_pointer<T>::value>> {
  using type = typename value_type<std::remove_pointer<T>>::type;
};

/** \ingroup type_trait
 * Helper function for accessing underlying type.
 */
template <typename T>
using value_type_t = typename value_type<T>::type;





#define STAN_ADD_REQUIRE_UNARY_VALUE(check_type, checker, doxygen_group) \
template <typename T> \
using require_##check_type##_vt = require_t<checker<value_type_t<std::decay_t<T>>>>; \
template <typename T> \
using require_not_##check_type##_vt = require_not_t<checker<value_type_t<std::decay_t<T>>>>;\
template <typename... Types> \
using require_all_##check_type##_vt = \
 require_all_t<checker<value_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_any_##check_type##_vt = \
 require_any_t<checker<value_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_all_not_##check_type##_vt \
    = require_all_not_t<checker<value_type_t<std::decay_t<Types>>>...>; \
template <typename... Types> \
using require_any_not_##check_type##_vt \
    = require_any_not_t<checker<value_type_t<std::decay_t<Types>>>...>; \


#define STAN_ADD_REQUIRE_BINARY_VALUE(check_type, checker) \
template <typename T, typename S> \
using require_##check_type##_vt = require_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>; \
template <typename T, typename S> \
using require_not_##check_type##_vt = require_not_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<S>>>>;\
template <typename T, typename... Types> \
using require_all_##check_type##_vt = \
 require_all_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>; \
template <typename T, typename... Types> \
using require_any_##check_type##_vt = \
 require_any_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>; \
template <typename T, typename... Types> \
using require_all_not_##check_type##_vt \
    = require_all_not_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>; \
template <typename T, typename... Types> \
using require_any_not_##check_type##_vt \
    = require_any_not_t<checker<value_type_t<std::decay_t<T>>, value_type_t<std::decay_t<Types>>>...>; \


}  // namespace stan

#endif
