#ifndef STAN_MATH_REV_META_IS_VARI_HPP
#define STAN_MATH_REV_META_IS_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a
 * vari.
 */
template <typename T>
struct is_vari<
    T, std::enable_if_t<std::is_base_of<
           math::vari_base, std::remove_pointer_t<std::decay_t<T>>>::value>>
    : std::true_type {};

namespace internal {
template <typename T, typename = void>
struct get_vari_value {
  using type = value_type_t<T>;
};

// until we figure out how to get inner type for vari_value
template <typename T>
struct get_vari_value<T, std::enable_if_t<is_vari<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};
}  // namespace internal

template <typename T>
using get_vari_t = typename internal::get_vari_value<T>::type;

template <template <class...> class TypeCheck, class... Check>
using require_vari_vt = require_t<
    container_type_check_base<is_vari, get_vari_t, TypeCheck, Check...>>;

template <template <class...> class TypeCheck, class... Check>
using require_all_vari_vt = require_all_t<
    container_type_check_base<is_vari, get_vari_t, TypeCheck, Check>...>;

}  // namespace stan
#endif
