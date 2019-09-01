#ifndef STAN_MATH_PRIM_MAT_META_GET_HPP
#define STAN_MATH_PRIM_MAT_META_GET_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>

namespace stan {

template <typename T, require_eigen<T>...>
inline auto&& get(T&& m, size_t n) {
  return std::forward<decltype(m.coeff(static_cast<int>(n)))>(m.coeff(static_cast<int>(n)));
}

}  // namespace stan
#endif
