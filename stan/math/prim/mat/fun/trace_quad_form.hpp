#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Compute trace(B^T A B).
 *
 * @tparam RA number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A matrix
 * @param B matrix
 * @return The trace of B^T A B
 * @throw std::domain_error if A is not square
 * @throw std::domain_error if A cannot be multiplied by B
 */
template <int RA, int CA, int RB, int CB>
inline double trace_quad_form(const Eigen::Matrix<double, RA, CA> &A,
                              const Eigen::Matrix<double, RB, CB> &B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  return (B.transpose() * A * B).trace();
}

}  // namespace math
}  // namespace stan

#endif
