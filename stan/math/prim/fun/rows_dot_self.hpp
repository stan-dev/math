#ifndef STAN_MATH_PRIM_FUN_ROWS_DOT_SELF_HPP
#define STAN_MATH_PRIM_FUN_ROWS_DOT_SELF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of each row of a matrix with itself.
 *
 * @tparam T type of the matrix (must be derived from \c
 * Eigen::MatrixBase)
 *
 * @param x matrix
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_st_var<T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, T::RowsAtCompileTime, 1> rows_dot_self(
    const T& x) {
  return x.rowwise().squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
