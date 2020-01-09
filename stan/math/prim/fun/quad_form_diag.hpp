#ifndef STAN_MATH_PRIM_FUN_QUAD_FORM_DIAG_HPP
#define STAN_MATH_PRIM_FUN_QUAD_FORM_DIAG_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, int R, int C>
inline Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, Eigen::Dynamic>
quad_form_diag(const Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>& mat,
               const Eigen::Matrix<T2, R, C>& vec) {
  check_vector("quad_form_diag", "vec", vec);
  check_square("quad_form_diag", "mat", mat);
  check_size_match("quad_form_diag", "rows of mat", mat.rows(), "size of vec",
                   vec.size());
  return vec.asDiagonal() * mat * vec.asDiagonal();
}

}  // namespace math
}  // namespace stan

#endif
