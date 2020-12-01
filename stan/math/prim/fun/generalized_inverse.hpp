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
  if (M.size() == 0) {
    return {};
  }

  constexpr int n = EigMat::RowsAtCompileTime;
  constexpr int m = EigMat::ColsAtCompileTime;

  if (n == m) {
      return M.inverse();
  }

  bool transp(false);
  constexpr int mn;

  if (n < m) {
    transp = true;
    mn = n;
    Eigen::Matrix<T, mn, mn> A = tcrossprod(G);
  } else {
     mn = m;
    Eigen::Matrix<T, mn, mn> A = = crossprod(G);
  }
  
   A = add_diag(A, rep_vector(1e-10, mn));

   Eigen::Matrix<T, mn, mn> L = cholesky_decompose(A);
   Eigen::Matrix<T, mn, mn> M = chol2inv(L);

   if (transp){
      return transpose(G) * tcrossprod(L * M);
   } else {
      return tcrossprod(L * M) * transpose(G);
   }
}

}  // namespace math
}  // namespace stan

#endif
