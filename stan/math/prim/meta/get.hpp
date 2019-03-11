#ifndef STAN_MATH_PRIM_META_GET_HPP
#define STAN_MATH_PRIM_META_GET_HPP

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {

template <typename T>
inline T get(const T& x, size_t n) {
  return x;
}

}  // namespace stan

namespace stan {

template <typename T>
inline T get(const std::vector<T>& x, size_t n) {
  return x[n];
}

}  // namespace stan

namespace stan {

template <typename T, int R, int C>
inline T get(const Eigen::Matrix<T, R, C>& m, size_t n) {
  return m(static_cast<int>(n));
}

}  // namespace stan
#endif
