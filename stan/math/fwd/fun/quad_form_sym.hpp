#ifndef STAN_MATH_FWD_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_FWD_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/fun/dot_product.hpp>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam TA type of elements in the symmetric matrix
 * @tparam RA number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam TB type of elements in the second matrix
 * @tparam RB number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam CB number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size CB.
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename TA, int RA, int CA, typename TB, int RB, int CB,
          require_any_fvar_t<TA, TB>...>
inline Eigen::Matrix<return_type_t<TA, TB>, CB, CB> quad_form_sym(
    const Eigen::Matrix<TA, RA, CA>& A, const Eigen::Matrix<TB, RB, CB>& B) {
  using T = return_type_t<TA, TB>;
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  Eigen::Matrix<T, CB, CB> ret(multiply(transpose(B), multiply(A, B)));
  return T(0.5) * (ret + transpose(ret));
}

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * @tparam TA type of elements in the symmetric matrix
 * @tparam RA number of rows in the symmetric matrix, can be Eigen::Dynamic
 * @tparam CA number of columns in the symmetric matrix, can be Eigen::Dynamic
 * @tparam TB type of elements in the vector
 * @tparam RB number of rows in the vector, can be Eigen::Dynamic
 *
 * @param A symmetric matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename TA, int RA, int CA, typename TB, int RB,
          require_any_fvar_t<TA, TB>...>
inline return_type_t<TA, TB> quad_form_sym(const Eigen::Matrix<TA, RA, CA>& A,
                                           const Eigen::Matrix<TB, RB, 1>& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  return dot_product(B, multiply(A, B));
}

}  // namespace math
}  // namespace stan

#endif
