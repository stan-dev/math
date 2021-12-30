#ifndef STAN_MATH_PRIM_FUN_NORM2_HPP
#define STAN_MATH_PRIM_FUN_NORM2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 * @return L2 norm of v.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_var, T>* = nullptr>
inline value_type_t<T> norm2(const T& v) {
  return v.template lpNorm<2>();
}

template <typename T>
inline T norm2(const std::vector<T>& x) {
  Eigen::Map<const Eigen::Matrix<T, -1, 1>> v(x.data(), x.size());
  return norm2(v);
}

}  // namespace math
}  // namespace stan

#endif
