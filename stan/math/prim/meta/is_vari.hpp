#ifndef STAN_MATH_PRIM_META_IS_VARI_HPP
#define STAN_MATH_PRIM_META_IS_VARI_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a
 * vari.
 */
template <typename T, typename = void>
struct is_vari : std::false_type {};

template <typename T>
struct value_type<T, require_t<is_vari<T>>> {
  using type = typename std::decay_t<T>::value_type;
};

}  // namespace stan
#endif
