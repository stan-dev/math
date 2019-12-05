#ifndef STAN_MATH_PRIM_MAT_META_GET_HPP
#define STAN_MATH_PRIM_MAT_META_GET_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {

template <typename T, typename = require_eigen_t<T>, typename = void>
inline auto get(const T& m, size_t n) {
  return m(static_cast<int>(n));
}

}  // namespace stan
#endif
