#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/adjoint_of.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>

namespace stan {
namespace math {

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param A first expression
 * @param B second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto multiply(T_a&& A, T_b&& B) {
  check_size_match("multiply ((OpenCL))", "A.cols()", A.cols(), "B.rows()",
                   B.rows());
  arena_t<T_a> a_arena = std::forward<T_a>(A);
  arena_t<T_b> b_arena = std::forward<T_b>(B);

  return make_callback_var(
      value_of(a_arena) * value_of(b_arena),
      [a_arena, b_arena](vari_value<matrix_cl<double>>& res) mutable {
        if (!is_constant<T_a>::value) {
          adjoint_of(a_arena) += res.adj() * transpose(value_of(b_arena));
        }
        if (!is_constant<T_b>::value) {
          adjoint_of(b_arena) += transpose(value_of(a_arena)) * res.adj();
        }
      });
}

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param A first expression
 * @param B second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto operator*(const T_a& A, const T_b& B) {
  return multiply(A, B);
}

/**
 * Return matrix multiplied by a scalar.
 *
 * @tparam T1 type of the scalar
 * @tparam T2 type of the matrix or expression
 *
 * @param A scalar
 * @param B matrix
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_stan_scalar_t<T1>* = nullptr,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T2>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline auto multiply(const T1& A, T2&& B) {
  arena_t<T2> b_arena = std::forward<T2>(B);

  return make_callback_var(
      value_of(A) * value_of(b_arena),
      [A, b_arena](vari_value<matrix_cl<double>>& res) mutable {
        adjoint_results(A, b_arena)
            += expressions(elt_multiply(res.adj(), value_of(b_arena)),
                           value_of(A) * res.adj());
      });
}

/**
 * Return matrix multiplied by a scalar.
 *
 * @tparam T1 type of the matrix or expression
 * @tparam T2 type of the scalar
 *
 * @param A matrix
 * @param B scalar
 * @return product of matrix and scalar
 */
template <typename T1, typename T2, require_stan_scalar_t<T2>* = nullptr,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T1>* = nullptr,
          require_any_var_t<T1, T2>* = nullptr>
inline auto multiply(const T1& A, const T2& B) {
  return multiply(B, A);
}

}  // namespace math
}  // namespace stan

#endif
#endif
