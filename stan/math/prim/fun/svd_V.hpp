#ifndef STAN_MATH_PRIM_FUN_SVD_V_HPP
#define STAN_MATH_PRIM_FUN_SVD_V_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Given input matrix m, return matrix V where `m = UDV^{T}`
 *
 * @tparam EigMat type of the matrix
 * @param m MxN input matrix
 * @return Orthogonal matrix V
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic, Eigen::Dynamic> svd_V(
    const EigMat& m) {
  check_nonzero_size("svd_V", "m", m);

  return Eigen::JacobiSVD<Eigen::Matrix<value_type_t<EigMat>, Eigen::Dynamic,
                                        Eigen::Dynamic> >(m,
                                                          Eigen::ComputeThinV)
      .matrixV();
}

}  // namespace math
}  // namespace stan

#endif
