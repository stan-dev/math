#ifndef STAN_MATH_REV_FUN_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun.hpp>
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
  if (!is_constant<T2>::value && !is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    auto arena_A_val = to_arena(arena_A.val());
    auto arena_B_val = to_arena(arena_B.val());
    using return_t
        = return_var_matrix_t<decltype(arena_A_val * arena_B_val), T1, T2>;
    arena_t<return_t> res = arena_A_val * arena_B_val;

    reverse_pass_callback(
        [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
          if (is_var_matrix<T1>::value || is_var_matrix<T2>::value) {
            arena_A.adj() += res.adj_op() * arena_B_val.transpose();
            arena_B.adj() += arena_A_val.transpose() * res.adj_op();
          } else {
            auto res_adj = res.adj().eval();
            arena_A.adj() += res_adj * arena_B_val.transpose();
            arena_B.adj() += arena_A_val.transpose() * res_adj;
          }
        });
    return return_t(res);
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<double, T1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    using return_t
        = return_var_matrix_t<decltype(arena_A * value_of(B).eval()), T1, T2>;
    arena_t<return_t> res = arena_A * arena_B.val_op();
    reverse_pass_callback([arena_B, arena_A, res]() mutable {
      arena_B.adj() += arena_A.transpose() * res.adj_op();
    });
    return return_t(res);
  } else {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<double, T2>> arena_B = value_of(B);
    using return_t
        = return_var_matrix_t<decltype(value_of(arena_A).eval() * arena_B), T1,
                              T2>;
    arena_t<return_t> res = arena_A.val_op() * arena_B;
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      arena_A.adj() += res.adj_op() * arena_B.transpose();
    });

    return return_t(res);
  }
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
  if (!is_constant<T2>::value && !is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<promote_scalar_t<double, T1>> arena_A_val = value_of(arena_A);
    arena_t<promote_scalar_t<double, T2>> arena_B_val = value_of(arena_B);
    var res = arena_A_val.dot(arena_B_val);

    reverse_pass_callback(
        [arena_A, arena_B, arena_A_val, arena_B_val, res]() mutable {
          auto res_adj = res.adj();
          arena_A.adj().array() += res_adj * arena_B_val.transpose().array();
          arena_B.adj().array() += arena_A_val.transpose().array() * res_adj;
        });
    return res;
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<promote_scalar_t<double, T1>> arena_A_val = value_of(A);
    var res = arena_A_val.dot(value_of(arena_B));
    reverse_pass_callback([arena_B, arena_A_val, res]() mutable {
      arena_B.adj().array() += arena_A_val.transpose().array() * res.adj();
    });
    return res;
  } else {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<double, T2>> arena_B_val = value_of(B);
    var res = value_of(arena_A).dot(arena_B_val);
    reverse_pass_callback([arena_A, arena_B_val, res]() mutable {
      arena_A.adj().array() += res.adj() * arena_B_val.transpose().array();
    });
    return res;
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
          require_return_type_t<is_var, T1, T2>* = nullptr,
          require_not_row_and_col_vector_t<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  if (!is_constant<T2>::value && !is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    using return_t = return_var_matrix_t<T2, T1, T2>;
    arena_t<return_t> res = arena_A.val() * arena_B.val().array();
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      const auto a_val = arena_A.val();
      for (Eigen::Index j = 0; j < res.cols(); ++j) {
        for (Eigen::Index i = 0; i < res.rows(); ++i) {
          const auto res_adj = res.adj().coeffRef(i, j);
          arena_A.adj() += res_adj * arena_B.val().coeff(i, j);
          arena_B.adj().coeffRef(i, j) += a_val * res_adj;
        }
      }
    });
    return return_t(res);
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<double, T1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    using return_t = return_var_matrix_t<T2, T1, T2>;
    arena_t<return_t> res = arena_A * arena_B.val().array();
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      arena_B.adj().array() += arena_A * res.adj().array();
    });
    return return_t(res);
  } else {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<double, T2>> arena_B = value_of(B);
    using return_t = return_var_matrix_t<T2, T1, T2>;
    arena_t<return_t> res = arena_A.val() * arena_B.array();
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      arena_A.adj() += (res.adj().array() * arena_B.array()).sum();
    });
    return return_t(res);
  }
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
inline auto multiply(const T1& A, const T2& B) {
  return multiply(B, A);
}

}  // namespace math
}  // namespace stan
#endif
