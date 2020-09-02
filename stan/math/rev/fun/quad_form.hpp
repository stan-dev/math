#ifndef STAN_MATH_REV_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam T1 type of the first (square) matrix
 * @tparam T2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline Eigen::Matrix<var, T2::ColsAtCompileTime, T2::ColsAtCompileTime>
quad_form_impl(const T1& A, const T2& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  using A_ref_t = ref_type_t<T1>;
  using B_ref_t = ref_type_t<T2>;

  A_ref_t A_ref = A;
  B_ref_t B_ref = B;

  check_not_nan("multiply", "A", A_ref);
  check_not_nan("multiply", "B", B_ref);

  arena_matrix<promote_scalar_t<double, T1>> arena_A_val;
  arena_matrix<promote_scalar_t<double, T2>> arena_B_val = value_of(B_ref);

  arena_matrix<promote_scalar_t<var, T1>> arena_A;
  arena_matrix<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
    arena_A_val = value_of(A_ref);
  }

  arena_matrix<Eigen::Matrix<var, T2::ColsAtCompileTime, T2::ColsAtCompileTime>>
      res;

  if(is_constant<T2>::value) {
    res = arena_B_val.transpose() * value_of(A_ref) * arena_B_val;
  } else {
    res = arena_B_val.transpose() * arena_A_val * arena_B_val;
  }

  reverse_pass_callback(
      [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
        auto C_adj = res.adj().eval();
	auto C_adj_B_t = (C_adj * arena_B_val.transpose()).eval();
	
        if (!is_constant<T1>::value)
          arena_A.adj() += arena_B_val * C_adj_B_t;

        if (!is_constant<T2>::value)
          arena_B.adj() += arena_A_val * C_adj_B_t.transpose() +
	    arena_A_val.transpose() * arena_B_val * C_adj;
      });

  return res;
}

}  // namespace internal

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam EigMat1 type of the first (square) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline auto quad_form(const EigMat1& A, const EigMat2& B) {
  return internal::quad_form_impl(A, B);
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam EigMat type of the matrix
 * @tparam ColVec type of the vector
 *
 * @param A square matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat, typename ColVec,
	  require_eigen_t<EigMat>* = nullptr,
          require_eigen_col_vector_t<ColVec>* = nullptr,
          require_any_vt_var<EigMat, ColVec>* = nullptr>
inline var quad_form(const EigMat& A, const ColVec& B) {
  return internal::quad_form_impl(A, B)(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
