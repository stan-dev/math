#ifndef STAN_MATH_PRIM_SCAL_META_GET_HPP
#define STAN_MATH_PRIM_SCAL_META_GET_HPP

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace stan {

template <typename T,
          typename = std::enable_if_t<std::is_fundamental<T>::value>>
inline T get(const T& x, size_t n) {
  return x;
}

}  // namespace stan
#endif
