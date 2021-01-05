#ifndef STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/multi_expression.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/multiply.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
auto elt_multiply(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  using inner_ret_type = decltype(value_of(m1).cwiseProduct(value_of(m2)));
  using ret_type = return_var_matrix_t<inner_ret_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.val().cwiseProduct(arena_m2.val()));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      auto ret_adj = ret.adj().array();
      eigen_results(arena_m1.adj().array(), arena_m2.adj().array())
          += eigen_expressions(arena_m2.val().array() * ret_adj,
                               arena_m1.val().array() * ret_adj);
    });
    return ret_type(ret);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    arena_t<ret_type> ret(arena_m1.val().cwiseProduct(arena_m2));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj().array() += arena_m2.array() * ret.adj().array();
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.cwiseProduct(arena_m2.val()));
    reverse_pass_callback([ret, arena_m2, arena_m1]() mutable {
      arena_m2.adj().array() += arena_m1.array() * ret.adj().array();
    });
    return ret_type(ret);
  }
}

/**
 * Return specified matrix multiplied by specified scalar where at least one
 * input has a scalar type of a `var_value`.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 *
 * @param A scalar
 * @param B matrix
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_not_matrix_t<T1>* = nullptr,
          require_matrix_t<T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto elt_multiply(const T1& A, const T2& B) {
  return multiply(A, B);
}

/**
 * Return specified matrix multiplied by specified scalar where at least one
 * input has a scalar type of a `var_value`.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 *
 * @param A matrix
 * @param B scalar
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_matrix_t<T1>* = nullptr,
          require_not_matrix_t<T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto elt_multiply(const T1& A, const T2& B) {
  return multiply(B, A);
}

}  // namespace math
}  // namespace stan

#endif
