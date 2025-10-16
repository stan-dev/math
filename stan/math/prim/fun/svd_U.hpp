#ifndef STAN_MATH_PRIM_FUN_SVD_U_HPP
#define STAN_MATH_PRIM_FUN_SVD_U_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Given input matrix m, return matrix U where `m = UDV^{T}`
 *
 * @tparam EigMat type of the matrix
 * @param m MxN input matrix
 * @return Orthogonal matrix U
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> svd_U(
    const EigMat& m) {
  using MatType
      = Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic>;
  if (unlikely(m.size() == 0)) {
    return MatType(0, 0);
  }
  return Eigen::JacobiSVD<MatType>(m, Eigen::ComputeThinU).matrixU();
}

}  // namespace math
}  // namespace stan

#endif
