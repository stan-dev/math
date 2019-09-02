#ifndef STAN_MATH_PRIM_SCAL_META_BOOL_CONSTANT_HPP
#define STAN_MATH_PRIM_SCAL_META_BOOL_CONSTANT_HPP

#include <type_traits>

namespace stan {
/**
 * Alias for structs used for wraps a static constant of bool.
 * @tparam B On true, inherits std::true_type, false is std::false_type
 */
template <bool B>
using bool_constant = std::integral_constant<bool, B>;
}  // namespace stan

#endif
