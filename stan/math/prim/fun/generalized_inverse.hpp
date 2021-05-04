#ifndef STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inverse.hpp>

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
  using value_t = value_type_t<EigMat>;

  if (G.size() == 0)
    return {};

  if (G.rows() == G.cols())
    return inverse(G);

  const auto& G_ref = to_ref(G);

  if (G.rows() < G.cols()) {
    return (G_ref * G_ref.transpose()).ldlt().solve(G_ref).transpose();
  } else {
    return (G_ref.transpose() * G_ref).ldlt().solve(G_ref.transpose());
  }
}

}  // namespace math
}  // namespace stan

#endif
