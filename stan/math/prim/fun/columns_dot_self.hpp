#ifndef STAN_MATH_PRIM_FUN_COLUMNS_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_COLUMNS_DOT_SELF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each column of a matrix with itself.
 *
 * @tparam T type of the matrix (must be derived from \c Eigen::MatrixBase)
 *
 * @param x Matrix.
 * @return Row vector containing the dot product of each column of the matrix
 * with itself.
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_var, T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, 1, T::ColsAtCompileTime> columns_dot_self(
    const T& x) {
  return x.colwise().squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
