#ifndef STAN_MATH_PRIM_SCAL_IS_STRING_CONVERTIBLE_HPP
#define STAN_MATH_PRIM_SCAL_IS_STRING_CONVERTIBLE_HPP

#include <type_traits>
#include <string>

namespace stan {

/**
 * Deduces whether type is convertible to string
 * @tparam T type to check
 */
template <typename T>
using is_string_convertible = std::is_convertible<T, std::string>;
}  // namespace stan

#endif
