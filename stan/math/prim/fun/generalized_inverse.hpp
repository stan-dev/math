#ifndef STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inverse.hpp>
#include <stan/math/prim/fun/generalized_inverse.hpp>

namespace stan {
namespace math {

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

  if (G_ref.rows() == G_ref.cols()) {
    Eigen::CompleteOrthogonalDecomposition<
        Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                      EigMat::ColsAtCompileTime>>
        complete_ortho_decomp_G = G_ref.completeOrthogonalDecomposition();
    if (!(complete_ortho_decomp_G.rank() < G_ref.rows()))
      return inverse(G_ref);
    else
      return complete_ortho_decomp_G.pseudoInverse();
  }

  if (G_ref.rows() < G_ref.cols()) {
    return (G_ref * G_ref.transpose()).ldlt().solve(G_ref).transpose();
  } else {
    return (G_ref.transpose() * G_ref).ldlt().solve(G_ref.transpose());
  }
}

}  // namespace math
}  // namespace stan

#endif
