#ifndef STAN_MATH_REV_META_IS_VAR_HPP
#define STAN_MATH_REV_META_IS_VAR_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier is a var.
 */
template <typename T>
struct is_var<T,
              std::enable_if_t<std::is_same<math::var, std::decay_t<T>>::value>>
    : std::true_type {};

namespace internal {

/** \ingroup type_trait
 * This underlying implementation is used when the type is not an std vector.
 */
template <typename T>
struct is_var_value_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implementation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_var_value_impl<math::var_value<Args...>> : std::true_type {};

}  // namespace internal

template <typename T>
struct is_var_value<
    T, std::enable_if_t<internal::is_var_value_impl<std::decay_t<T>>::value>> {
};

namespace internal {
template <typename T, typename = void>
struct get_var_value {
  using type = value_type_t<T>;
};

// until we figure out how to get inner type for vari_value
template <typename T>
struct get_var_value<T, std::enable_if_t<is_var<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};
template <typename T, typename = void>
struct get_var_vari_value {
  using type = value_type_t<T>;
};

// until we figure out how to get inner type for vari_value
template <typename T>
struct get_var_vari_value<T, std::enable_if_t<is_var_value<T>::value>> {
  using type = typename std::decay_t<T>::vari_type;
};

template <typename T>
using get_var_vari_value_t = typename get_var_vari_value<std::decay_t<T>>::type;

}  // namespace internal


}  // namespace stan
#endif
