#ifndef STAN_MATH_PRIM_SCAL_META_GET_HPP
#define STAN_MATH_PRIM_SCAL_META_GET_HPP

#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <cmath>
#include <cstddef>

namespace stan {

template <typename T, require_stan_scalar<T>...>
inline auto&& get(T&& x, size_t n) {
  return std::forward<T>(x);
}

}  // namespace stan
#endif
