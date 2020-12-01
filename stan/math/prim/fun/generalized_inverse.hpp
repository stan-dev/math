#ifndef STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/chol2inv.hpp>
#include <stan/math/prim/fun/add_diag.hpp>
#include <stan/math/prim/fun/rep_vector.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>
#include <stan/math/prim/fun/crossprod.hpp>

namespace stan {
namespace math {

/**
 * Returns the Moore-Penrose generalized inverse of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam n number of rows, can be Eigen::Dynamic
 * @tparam m number of columns, can be Eigen::Dynamic
 *
 * @param M specified matrix
 * @return Generalized inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 */
template <typename EigMat, require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
generalized_inverse(const EigMat& M) {
  using value_t = value_type_t<EigMat>;
  if (M.size() == 0) {
    return {};
  }
  const auto n = M.rows();
  const auto m = M.cols();

  if (n == m) {
      return M.inverse();
  } else if (n < m) {
    Eigen::Matrix<value_t, -1, -1> A = tcrossprod(M);
    A.diagonal().array() += Eigen::Array<double, -1, 1>::Constant(n, 1e-10);
    Eigen::Matrix<value_t, -1, -1> L = cholesky_decompose(A);
    Eigen::Matrix<value_t, -1, -1> M = chol2inv(L);
    return transpose(M) * tcrossprod(L * M);
  } else {
    Eigen::Matrix<value_t, -1, -1> A = crossprod(M);
    A.diagonal().array() += Eigen::Array<double, -1, 1>::Constant(m, 1e-10);
    Eigen::Matrix<value_t, -1, -1> L = cholesky_decompose(A);
    Eigen::Matrix<value_t, -1, -1> M = chol2inv(L);
    return tcrossprod(L * M) * transpose(M);
  }
}

}  // namespace math
}  // namespace stan

#endif
