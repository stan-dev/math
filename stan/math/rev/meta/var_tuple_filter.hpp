#ifndef STAN_MATH_REV_META_VAR_TO_VARI_FILTER_HPP
#define STAN_MATH_REV_META_VAR_TO_VARI_FILTER_HPP

#include <stan/math/rev/meta/is_var.hpp>
namespace stan {
// NOTE: This should probably be in Stan namespace not math
namespace math {

namespace internal {
template <typename T>
using container_var_vari_value_t
    = std::conditional_t<is_container<std::decay_t<T>>::value, vari**, vari*>;

template <typename T>
using contains_var_value = is_var<scalar_type_t<T>>;
}  // namespace internal

/**
 * Given a paramter pack of types, return a tuple with `vari*` for var types and
 * `vari**` for containers of var types. non-var types will have an element
 * equal to `nullptr`.
 * tparam Ts Any set of types.
 */
template <typename... Ts>
using var_to_vari_filter_t = decltype(std::tuple_cat(std::declval<
    std::conditional_t<internal::contains_var_value<Ts>::value,
      std::tuple<internal::container_var_vari_value_t<Ts>>,
      std::tuple<std::nullptr_t>>>()...));

}  // namespace math
}  // namespace stan
#endif
