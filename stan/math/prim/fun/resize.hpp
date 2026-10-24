#ifndef STAN_MATH_PRIM_FUN_RESIZE_HPP
#define STAN_MATH_PRIM_FUN_RESIZE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename T, int R, int C>
inline void resize(Eigen::Matrix<T, R, C>& x, const std::vector<int>& dims,
                   int pos) {
  if constexpr (R == 1 || C == 1) {
    if (pos + 1 == static_cast<int>(dims.size())) {
      if constexpr (R == 1 && C != 1) {
        x.resize(1, dims[pos]);
      } else {
        x.resize(dims[pos], 1);
      }
      return;
    }
  }
  x.resize(dims[pos], dims[pos + 1]);
}

template <typename T>
inline void resize(T /*x*/, const std::vector<int>& /*dims*/,
                   int /*pos*/ = 0) {  // no-op
}

template <typename T>
inline void resize(std::vector<T>& x, const std::vector<int>& dims,
                   int pos = 0) {
  x.resize(dims[pos]);
  ++pos;

  if (pos >= static_cast<int>(dims.size())) {
    return;  // skips lowest loop to scalar
  }

  for (size_t i = 0; i < x.size(); ++i) {
    resize(x[i], dims, pos);
  }
}
}  // namespace internal

/**
 * Recursively resize the specified vector of vectors,
 * which must bottom out at scalar values, Eigen vectors
 * or Eigen matrices.
 *
 * @tparam T type of object being resized
 * @param x Array-like object to resize.
 * @param dims New dimensions.
 */
template <typename T>
inline void resize(T& x, std::vector<int> dims) {
  internal::resize(x, dims, 0U);
}

}  // namespace math
}  // namespace stan
#endif
