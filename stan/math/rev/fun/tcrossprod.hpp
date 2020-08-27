#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam T Type of the matrix (must be derived from \c Eigen::MatrixBase)
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline Eigen::Matrix<var, T::RowsAtCompileTime, T::RowsAtCompileTime>
tcrossprod(const T& M) {
  if (M.rows() == 0) {
    return {};
  }
  // if (M.rows() == 1)
  //   return M * M.transpose();

  // WAS JUST THIS
  // matrix_v result(M.rows(), M.rows());
  // return result.setZero().selfadjointView<Eigen::Upper>().rankUpdate(M);

  Eigen::Matrix<var, T::RowsAtCompileTime, T::RowsAtCompileTime> MMt(M.rows(),
                                                                     M.rows());

  vari** vs
      = reinterpret_cast<vari**>(ChainableStack::instance_->memalloc_.alloc(
          (M.rows() * M.cols()) * sizeof(vari*)));
  int pos = 0;
  for (int m = 0; m < M.rows(); ++m) {
    for (int n = 0; n < M.cols(); ++n) {
      vs[pos++] = M(m, n).vi_;
    }
  }
  for (int m = 0; m < M.rows(); ++m) {
    MMt.coeffRef(m, m)
        = var(new internal::dot_self_vari(vs + m * M.cols(), M.cols()));
  }
  for (int m = 0; m < M.rows(); ++m) {
    for (int n = 0; n < m; ++n) {
      MMt.coeffRef(n, m) = MMt.coeffRef(m, n) = dot_product(M.row(m), M.row(n));
    }
  }
  return MMt;
}

}  // namespace math
}  // namespace stan
#endif
