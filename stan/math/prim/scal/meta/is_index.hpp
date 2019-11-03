#ifndef STAN_MATH_PRIM_SCAL_META_IS_INDEX_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_INDEX_HPP

#include <type_traits>

namespace stan {

/**
 * Deduces whether type is non-floating point arithmetic
 * @tparam T type to check
 */
template <typename T>
using is_index = bool_constant<!std::is_floating_point<std::decay_t<T>>::value
                               && std::is_arithmetic<std::decay_t<T>>::value>;

}

#endif
