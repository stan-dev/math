#ifndef STAN_MATH_PRIM_FUN_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_DOT_SELF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T dot_self(const T& x) {
  return x * x;
}

inline double dot_self(const std::vector<double>& x) {
  double sum = 0.0;
  for (double i : x) {
    sum += i * i;
  }
  return sum;
}

/**
 * Returns squared norm of a vector or matrix. For vectors that equals the dot
 * product of the specified vector with itself.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_var, T>* = nullptr>
inline value_type_t<T> dot_self(const T& v) {
  return v.squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
