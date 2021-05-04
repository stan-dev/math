#ifndef STAN_MATH_PRIM_FUN_TCROSSPROD_HPP
#define STAN_MATH_PRIM_FUN_TCROSSPROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline Eigen::Matrix<value_type_t<T>, T::RowsAtCompileTime,
                     T::RowsAtCompileTime>
tcrossprod(const T& M) {
  if (M.rows() == 0) {
    return {};
  }
  const auto& M_ref = to_ref(M);
  if (M.rows() == 1) {
    return M_ref * M_ref.transpose();
  }
  Eigen::Matrix<value_type_t<T>, T::RowsAtCompileTime, T::RowsAtCompileTime>
      result(M.rows(), M.rows());
  return result.setZero().template selfadjointView<Eigen::Upper>().rankUpdate(
      M_ref);
}

}  // namespace math
}  // namespace stan

#endif
