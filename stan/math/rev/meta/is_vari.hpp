#ifndef STAN_MATH_REV_META_IS_VARI_HPP
#define STAN_MATH_REV_META_IS_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
/** \ingroup type_trait
 * Specialization for checking if value of T minus cv qualifier and pointer is a vari.
 */
 template <typename T>
 struct is_vari : bool_constant<std::is_base_of<math::vari_base, std::remove_pointer_t<std::decay_t<T>>>::value> {};

template <typename T>
using require_vari_t = require_t<is_vari<T>>;

template <typename... Types>
using require_all_vari_t = require_all_t<is_vari<Types>...>;

}  // namespace stan
#endif
