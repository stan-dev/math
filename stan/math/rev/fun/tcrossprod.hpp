#ifndef STAN_MATH_REV_FUN_TCROSSPROD_HPP
#define STAN_MATH_REV_FUN_TCROSSPROD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/Eigen_NumTraits.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of post-multiplying a matrix by its
 * own transpose.
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param M Matrix to multiply.
 * @return M times its transpose.
 */
template <int R, int C>
inline Eigen::Matrix<var, -1, -1> tcrossprod(
    const Eigen::Matrix<var, R, C>& M) {
  if (M.rows() == 0) {
    return {};
  }
  // if (M.rows() == 1)
  //   return M * M.transpose();

  // WAS JUST THIS
  // matrix_v result(M.rows(), M.rows());
  // return result.setZero().selfadjointView<Eigen::Upper>().rankUpdate(M);

  matrix_v MMt(M.rows(), M.rows());

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
    MMt(m, m) = var(new internal::dot_self_vari(vs + m * M.cols(), M.cols()));
  }
  for (int m = 0; m < M.rows(); ++m) {
    for (int n = 0; n < m; ++n) {
      MMt(m, n) = var(new internal::dot_product_vari<var, var>(
          vs + m * M.cols(), vs + n * M.cols(), M.cols()));
      MMt(n, m) = MMt(m, n);
    }
  }
  return MMt;
}

}  // namespace math
}  // namespace stan
#endif
