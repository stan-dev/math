#ifndef STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, int R1, int C1, int R2, int C2>
Eigen::Matrix<return_type_t<T1, T2>, R1, C1> diag_post_multiply(
    const Eigen::Matrix<T1, R1, C1>& m1, const Eigen::Matrix<T2, R2, C2>& m2) {
  check_vector("diag_post_multiply", "m2", m2);
  check_size_match("diag_post_multiply", "m2.size()", m2.size(), "m1.cols()",
                   m1.cols());
  return m1 * m2.asDiagonal();
}

}  // namespace math
}  // namespace stan
#endif
