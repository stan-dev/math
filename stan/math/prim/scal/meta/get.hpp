#ifndef STAN_MATH_PRIM_SCAL_META_GET_HPP
#define STAN_MATH_PRIM_SCAL_META_GET_HPP

#include <stan/math/prim/arr/meta/get.hpp>
#include <cmath>
#include <cstddef>

namespace stan {

template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, T>* = nullptr>
inline T get(const T& x, size_t n) {
  return x;
}

}  // namespace stan
#endif
