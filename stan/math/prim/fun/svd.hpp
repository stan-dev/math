#ifndef STAN_MATH_PRIM_FUN_SVD_HPP
#define STAN_MATH_PRIM_FUN_SVD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Given input matrix m, return the singular value decomposition (U,D,V)
 * such that `m = U*diag(D)*V^{T}`
 *
 * @tparam EigMat type of the matrix
 * @param m MxN input matrix
 * @return a tuple (U,D,V) where U is an orthogonal matrix, D a vector of
 * singular values (in decreasing order), and V an orthogonal matrix
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_st_var<EigMat>* = nullptr>
std::tuple<Eigen::Matrix<value_type_t<EigMat>, -1, -1>,
           Eigen::Matrix<base_type_t<EigMat>, -1, 1>,
           Eigen::Matrix<value_type_t<EigMat>, -1, -1>>
svd(const EigMat& m) {
  if (unlikely(m.size() == 0)) {
    return std::make_tuple(Eigen::Matrix<value_type_t<EigMat>, -1, -1>(0, 0),
                           Eigen::Matrix<base_type_t<EigMat>, -1, 1>(0, 1),
                           Eigen::Matrix<value_type_t<EigMat>, -1, -1>(0, 0));
  }

  Eigen::JacobiSVD<Eigen::Matrix<value_type_t<EigMat>, -1, -1>> svd(
      m, Eigen::ComputeThinU | Eigen::ComputeThinV);
  return std::make_tuple(std::move(svd.matrixU()),
                         std::move(svd.singularValues()),
                         std::move(svd.matrixV()));
}

}  // namespace math
}  // namespace stan

#endif
