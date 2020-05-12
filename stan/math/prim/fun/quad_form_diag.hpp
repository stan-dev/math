#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_DIAG_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_DIAG_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename EigMat, typename EigVec, require_eigen_t<EigMat>* = nullptr,
          require_eigen_vector_t<EigVec>* = nullptr>
inline Eigen::Matrix<return_type_t<EigMat, EigVec>, Eigen::Dynamic,
                     Eigen::Dynamic>
quad_form_diag(const EigMat& mat, const EigVec& vec) {
  check_square("quad_form_diag", "mat", mat);
  check_size_match("quad_form_diag", "rows of mat", mat.rows(), "size of vec",
                   vec.size());
  return vec.asDiagonal() * mat * vec.asDiagonal();
}

}  // namespace math
}  // namespace stan

#endif
