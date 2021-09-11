#ifndef STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/crossprod.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inverse.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>

namespace stan {
namespace math {

namespace internal {
  template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::ColsAtCompileTime,
                     EigMat::RowsAtCompileTime>
 cholesky_low_rank_decomposition(const EigMat& A) {
    if (A.size() == 0)
    return {};
    
    int n = A.rows();
    auto dA = A.diagonal();
    double tol = 1e-12;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> L = Eigen::MatrixXd::Zero(n, n);
    int r = 0;
  
    for (int i = 0; i < n; i++) {
      if (r == 0) {
        L.block(i, r, n - i, 1) = A.block(i, i, n - i, 1);
      } else {
        L.block(i, r , n - i, 1) = A.block(i, i, n - i, 1) - (L.block(i, 0, n - i, r) * L.block(i, 0, 1, r).transpose());
      }

      if (L(i, r) > tol ) {
        L(i, r) = sqrt(L(i, r)) ;
        if (i < n - 1) {
             L.block(i + 1, r, n - (i + 1), 1) = (L.block(i + 1, r, n - (i + 1), 1).array() / L(i, r)).matrix(); 
        }
      }  else {
        r -= 1;
      }
      r += 1;
    }
  return L.block(0, 0, n, r);
 }
} // namespace internal

/**
 * Returns the Moore-Penrose generalized inverse of the specified matrix.
 *
 * The method is based on the Cholesky computation of the transform as specified
 *  in
 *
 * <ul><li> Courrieu, Pierre. 2008.  Fast Computation of Moore-Penrose Inverse
 Matrices.
 * <i>arXiv</i> <b>0804.4809</b> </li></ul>
 *
 * @tparam EigMat type of the matrix (must be derived from `Eigen::MatrixBase`)
 *
 * @param G specified matrix
 * @return Generalized inverse of the matrix (an empty matrix if the specified
 * matrix has size zero).
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr,
          require_not_vt_var<EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::ColsAtCompileTime,
                     EigMat::RowsAtCompileTime>
generalized_inverse(const EigMat& G) {
  const auto& G_ref = to_ref(G);
  if (G_ref.size() == 0) 
    return {};
  const int n = std::min(G.rows(), G.cols());
  const bool transpose_bool = G.rows() == n ? true : false;

  auto A = transpose_bool ? tcrossprod(G) : crossprod(G);
  // note: L may not be square. So L * M ops don't cancel.
  auto L = internal::cholesky_low_rank_decomposition(A);
  auto M = inverse(crossprod(L));

  if (transpose_bool) {
    return G_ref.transpose() * L * M * M * L.transpose();
  } else {
    return L * M * M * L.transpose() * G_ref.transpose();
  }
}

}  // namespace math
}  // namespace stan

#endif
