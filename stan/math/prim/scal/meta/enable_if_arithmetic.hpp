#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_ARITHMETIC_HPP

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_arithmetic = std::enable_if_t<std::is_arithmetic<T>::value>;

}  // namespace stan
#endif
