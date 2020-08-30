#ifndef STAN_MATH_REV_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {

/**
 * Return the product of two matrices.
 *
 * @tparam T1 type of first matrix
 * @tparam T2 type of second matrix
 *
 * @param[in] A first matrix
 * @param[in] B second matrix
 * @return A * B
 */
template <typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline Eigen::Matrix<var,
		     T1::RowsAtCompileTime,
		     T2::ColsAtCompileTime>
multiply_impl(const T1& A, const T2& B) {
  check_multiplicable("multiply", "A", A, "B", B);

  using A_ref_t = ref_type_t<T1>;
  using B_ref_t = ref_type_t<T2>;
  
  A_ref_t A_ref = A;
  B_ref_t B_ref = B;

  check_not_nan("multiply", "A", A_ref);
  check_not_nan("multiply", "B", B_ref);

  arena_matrix<promote_scalar_t<double, T1>> arena_A_val = value_of(A_ref);
  arena_matrix<promote_scalar_t<double, T2>> arena_B_val = value_of(B_ref);

  arena_matrix<promote_scalar_t<var, T1>> arena_A;
  arena_matrix<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_ref;
    arena_B_val = value_of(B_ref);
  }

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
    arena_A_val = value_of(A_ref);
  }

  arena_matrix<Eigen::Matrix<var,
			     T1::RowsAtCompileTime,
			     T2::ColsAtCompileTime>> res;

  if(!is_constant<T1, T2>::value) {
    res = arena_A_val * arena_B_val;
  } else if(!is_constant<T1>::value) {
    res = value_of(A_ref) * arena_B_val;
  } else if(!is_constant<T2>::value) {
    res = arena_A_val * value_of(B_ref);
  }

  reverse_pass_callback([arena_A, arena_B,
			 arena_A_val, arena_B_val,
			 res]() mutable {
    auto res_adj = res.adj().eval();

    if (!is_constant<T1>::value)
      arena_A.adj() += res_adj * arena_B_val.transpose();

    if (!is_constant<T2>::value)
      arena_B.adj() += arena_A_val.transpose() * res_adj;
  });

  return res;
}

}

/**
 * Return the product of two matrices.
 *
 * This version does not handle row vector times column vector
 *
 * @tparam T1 type of first matrix
 * @tparam T2 type of second matrix
 *
 * @param[in] A first matrix
 * @param[in] B second matrix
 * @return A * B
 */
template <typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr,
	  require_not_eigen_row_and_col_t<T1, T2>* = nullptr>
inline Eigen::Matrix<var,
		     T1::RowsAtCompileTime,
		     T2::ColsAtCompileTime>
multiply(const T1& A, const T2& B) {
  return internal::multiply_impl(A, B);
}

/**
 * Return the product of a row vector times a column vector as a scalar
 *
 * @tparam T1 type of row vector
 * @tparam T2 type of column vector
 *
 * @param[in] A row vector
 * @param[in] B column vector
 * @return A * B as a scalar
 */
template <typename T1, typename T2,
          require_all_eigen_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr,
	  require_eigen_row_and_col_t<T1, T2>* = nullptr>
inline var multiply(const T1& A, const T2& B) {
  return internal::multiply_impl(A, B)(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
