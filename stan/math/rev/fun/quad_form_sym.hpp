#ifndef STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form_sym.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam EigMat1 type of the first (symmetric) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return the quadratic form as a matrix
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<var, EigMat2::ColsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
quad_form_sym_impl(const EigMat1& A, const EigMat2& B) {
  check_multiplicable("quad_form_sym", "A", A, "B", B);

  using A_ref_t = ref_type_t<EigMat1>;
  using B_ref_t = ref_type_t<EigMat2>;

  A_ref_t A_ref = A;
  B_ref_t B_ref = B;

  check_symmetric("quad_form_sym", "A", A_ref);

  check_not_nan("multiply", "A", A_ref);
  check_not_nan("multiply", "B", B_ref);

  arena_matrix<promote_scalar_t<double, EigMat1>> arena_A_val;
  arena_matrix<promote_scalar_t<double, EigMat2>> arena_B_val = value_of(B_ref);

  arena_matrix<promote_scalar_t<var, EigMat1>> arena_A;
  arena_matrix<promote_scalar_t<var, EigMat2>> arena_B;

  if (!is_constant<EigMat1>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<EigMat2>::value) {
    arena_B = B_ref;
    arena_A_val = value_of(A_ref);
  }

  arena_matrix<Eigen::MatrixXd> res_val;

  if (is_constant<EigMat2>::value) {
    res_val = arena_B_val.transpose() * value_of(A_ref) * arena_B_val;
  } else {
    res_val = arena_B_val.transpose() * arena_A_val * arena_B_val;
  }

  arena_matrix<Eigen::Matrix<var, EigMat2::ColsAtCompileTime,
                             EigMat2::ColsAtCompileTime>>
      res = 0.5 * (res_val + res_val.transpose());

  reverse_pass_callback(
      [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
        auto C_adj = res.adj().eval();
        auto C_adj_B_t = (C_adj * arena_B_val.transpose()).eval();

        if (!is_constant<EigMat1>::value)
          arena_A.adj() += arena_B_val * C_adj_B_t;

        if (!is_constant<EigMat2>::value)
          arena_B.adj() += arena_A_val * C_adj_B_t.transpose()
                           + arena_A_val.transpose() * arena_B_val * C_adj;
      });

  return res;
}
}  // namespace internal

/**
 * Return the quadratic form \f$ B^T A B \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam EigMat1 type of the first (symmetric) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A symmetric matrix
 * @param B second matrix
 * @return the quadratic form as a matrix
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline auto quad_form_sym(const EigMat1& A, const EigMat2& B) {
  return internal::quad_form_sym_impl(A, B);
}

/**
 * Return the quadratic form \f$ b^T A b \f$ of a symmetric matrix.
 *
 * Symmetry of the resulting matrix is guaranteed.
 *
 * @tparam EigMat type of the symmetric matrix A
 * @tparam ColVec type of the vector b
 *
 * @param A symmetric matrix
 * @param b column vector
 * @return the quadratic form as a scalar
 * @throws std::invalid_argument if A is not symmetric, or if A cannot be
 * multiplied by b
 */
template <typename EigMat, typename ColVec, require_eigen_t<EigMat>* = nullptr,
          require_eigen_col_vector_t<ColVec>* = nullptr,
          require_any_vt_var<EigMat, ColVec>* = nullptr>
inline var quad_form_sym(const EigMat& A, const ColVec& b) {
  return internal::quad_form_sym_impl(A, b)(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
