#ifndef STAN_MATH_REV_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim.hpp>
#include <type_traits>

namespace stan {
namespace math {

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
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
         require_return_type_t<is_var, T1, T2>* = nullptr,
         require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  check_multiplicable("multiply", "A", A, "B", B);

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  check_not_nan("multiply", "A", A_ref);
  check_not_nan("multiply", "B", B_ref);


  arena_t<promote_scalar_t<var, T1>> arena_A = to_arena_if<!is_constant<T1>::value>(A_ref);
  arena_t<promote_scalar_t<var, T2>> arena_B = to_arena_if<!is_constant<T2>::value>(B_ref);
  arena_t<promote_scalar_t<double, T1>> arena_A_val = to_arena_if<!is_constant<T2>::value>(value_of(A_ref));
  arena_t<promote_scalar_t<double, T2>> arena_B_val = to_arena_if<!is_constant<T1>::value>(value_of(B_ref));
  using return_t = promote_var_matrix_t<decltype(arena_A_val * arena_B_val), T1, T2>;
  arena_t<return_t> res;

  if (!is_constant<T1>::value) {
    res = value_of(A_ref) * arena_B_val;
  } else if (!is_constant<T2>::value) {
    res = arena_A_val * value_of(B_ref);
  } else {
    res = arena_A_val * arena_B_val;
  }

  reverse_pass_callback(
      [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
        auto res_adj = res.adj().eval();

        if (!is_constant<T1>::value)
          arena_A.adj() += res_adj * arena_B_val.transpose();

        if (!is_constant<T2>::value)
          arena_B.adj() += arena_A_val.transpose() * res_adj;
      });

  return return_t(res);
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
 template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
           require_return_type_t<is_var, T1, T2>* = nullptr,
           require_row_and_col_vector_t<T1, T2>* = nullptr>
inline var multiply(const T1& A, const T2& B) {
  check_multiplicable("multiply", "A", A, "B", B);

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  check_not_nan("multiply", "A", A_ref);
  check_not_nan("multiply", "B", B_ref);

  arena_t<promote_scalar_t<var, T1>> arena_A = to_arena_if<!is_constant<T1>::value>(A_ref);
  arena_t<promote_scalar_t<var, T2>> arena_B = to_arena_if<!is_constant<T2>::value>(B_ref);
  arena_t<promote_scalar_t<double, T1>> arena_A_val = to_arena_if<!is_constant<T2>::value>(value_of(A_ref));
  arena_t<promote_scalar_t<double, T2>> arena_B_val = to_arena_if<!is_constant<T1>::value>(value_of(B_ref));
  arena_t<var> res;

  if (!is_constant<T1>::value) {
    res = value_of(A_ref).dot(arena_B_val);
  } else if (!is_constant<T2>::value) {
    res = arena_A_val.dot(value_of(B_ref));
  } else {
    res = arena_A_val.dot(arena_B_val);
  }

  reverse_pass_callback(
      [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
        auto res_adj = res.adj();

        if (!is_constant<T1>::value)
          arena_A.adj().array() += res_adj * arena_B_val.transpose().array();

        if (!is_constant<T2>::value)
          arena_B.adj().array() += arena_A_val.transpose().array() * res_adj;
      });

  return res;
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
          require_return_type_t<is_var, T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  const auto& B_ref = to_ref(B);

  check_not_nan("multiply", "A", A);
  check_not_nan("multiply", "B", B_ref);


  arena_t<promote_scalar_t<var, T2>> arena_B = to_arena_if<!is_constant<T2>::value>(B_ref);
  arena_t<promote_scalar_t<double, T2>> arena_B_val = to_arena_if<!is_constant<T1>::value>(value_of(B_ref));
  using return_t = promote_var_matrix_t<T2, T1, T2>;
  arena_t<return_t> res;

  if (!is_constant<T1>::value) {
    res = value_of(A) * arena_B_val.array();
  } else if (!is_constant<T2>::value) {
    res = value_of(A) * value_of(B_ref).array();
  } else {
    res = value_of(A) * arena_B_val.array();
  }

  reverse_pass_callback(
      [A, arena_B, arena_B_val, res]() mutable {
        auto res_adj = res.adj().eval();

        if (!is_constant<T1>::value)
          forward_as<var>(A).adj() += (res_adj.array() * arena_B_val.array()).sum();

        if (!is_constant<T2>::value)
          arena_B.adj().array() += value_of(A) * res_adj.array();
      });

  return return_t(res);
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
           require_return_type_t<is_var, T1, T2>* = nullptr,
           require_not_row_and_col_vector_t<T1, T2>* = nullptr>
 inline auto multiply(const T1& A, const T2& B) {
   const auto& A_ref = to_ref(A);

   check_not_nan("multiply", "A", A_ref);
   check_not_nan("multiply", "B", B);


   arena_t<promote_scalar_t<var, T1>> arena_A = to_arena_if<!is_constant<T1>::value>(A_ref);
   arena_t<promote_scalar_t<double, T1>> arena_A_val = to_arena_if<!is_constant<T2>::value>(value_of(A_ref));
   using return_t = promote_var_matrix_t<T1, T1, T2>;
   arena_t<return_t> res;

   if (!is_constant<T1>::value) {
     res = value_of(A_ref).array() * value_of(B);
   } else if (!is_constant<T2>::value) {
     res = arena_A_val.array() * value_of(B);
   } else {
     res = arena_A_val.array() * value_of(B);
   }

   reverse_pass_callback(
       [arena_A, B, arena_A_val, res]() mutable {
         auto res_adj = res.adj().eval();

         if (!is_constant<T1>::value)
           arena_A.adj().array() += value_of(B) * res_adj.array();

         if (!is_constant<T2>::value)
           forward_as<var>(B).adj() += (res_adj.array() * arena_A_val.array()).sum();
       });

   return return_t(res);
  }

}  // namespace math
}  // namespace stan
#endif
