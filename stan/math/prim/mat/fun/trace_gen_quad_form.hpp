#ifndef STAN_MATH_PRIM_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/trace.hpp>
#include <stan/math/prim/mat/fun/multiply.hpp>
#include <stan/math/prim/mat/fun/transpose.hpp>
#include <exception>

namespace stan {
namespace math {

/**
 * Return the trace of D times the quadratic form of B and A.
 * That is, `trace_gen_quad_form(D, A, B) = trace(D * B' * A * B).`
 *
 * @tparam TD type of elements in the first matrix
 * @tparam TA type of elements in the second matrix
 * @tparam TB type of elements in the third matrix
 * @tparam RD number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam CD number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam RA number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the second matrix, can be Eigen::Dynamic
 * @tparam RB number of rows in the third matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the third matrix, can be Eigen::Dynamic
 *
 * @param D multiplier
 * @param A outside term in quadratic form
 * @param B inner term in quadratic form
 * @return trace(D * B' * A * B)
 * @throw std::domain_error if A or D is not square
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by D.
 */
template <typename TD, int RD, int CD, typename TA, int RA, int CA, typename TB,
          int RB, int CB, typename = require_any_not_var_t<TD, TA, TB>>
inline return_type_t<TD, TA, TB> trace_gen_quad_form(
    const Eigen::Matrix<TD, RD, CD> &D, const Eigen::Matrix<TA, RA, CA> &A,
    const Eigen::Matrix<TB, RB, CB> &B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);
  return trace(multiply(multiply(D, transpose(B)), multiply(A, B)));
}

/**
 * Return the trace of D times the quadratic form of B and A.
 * That is, `trace_gen_quad_form(D, A, B) = trace(D * B' * A * B).`
 * This is the double-only overload to allow Eigen's expression
 * templates to be used for efficiency.
 *
 * @tparam RD number of rows in the first matrix, can be Eigen::Dynamic
 * @tparam CD number of columns in the first matrix, can be Eigen::Dynamic
 * @tparam RA number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the second matrix, can be Eigen::Dynamic
 * @tparam TB type of elements in the third matrix
 * @tparam RB number of rows in the third matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the third matrix, can be Eigen::Dynamic
 *
 * @param D multiplier
 * @param A outside term in quadratic form
 * @param B inner term in quadratic form
 * @return trace(D * B' * A * B)
 * @throw std::domain_error if A or D is not square
 * @throw std::domain_error if A cannot be multiplied by B or B cannot
 * be multiplied by D.
 */
template <int RD, int CD, int RA, int CA, typename TB, int RB, int CB>
inline double trace_gen_quad_form(const Eigen::Matrix<double, RD, CD> &D,
                                  const Eigen::Matrix<double, RA, CA> &A,
                                  const Eigen::Matrix<double, RB, CB> &B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);
  return (D * B.transpose() * A * B).trace();
}

}  // namespace math
}  // namespace stan

#endif
