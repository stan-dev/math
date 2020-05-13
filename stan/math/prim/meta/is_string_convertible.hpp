#ifndef STAN_MATH_PRIM_META_IS_STRING_CONVERTIBLE_HPP
#define STAN_MATH_PRIM_META_IS_STRING_CONVERTIBLE_HPP

#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <string>

namespace stan {

/**
 * Deduces whether type is convertible to string
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
using is_string_convertible = std::is_convertible<T, std::string>;

STAN_ADD_REQUIRE_UNARY(string_convertible, is_string_convertible, require_std);
STAN_ADD_REQUIRE_UNARY_INNER(string_convertible, is_string_convertible,
                             require_std);

}  // namespace stan

#endif
