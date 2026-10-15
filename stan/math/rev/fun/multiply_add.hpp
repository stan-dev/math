#ifndef STAN_MATH_REV_FUN_MULTIPLY_ADD_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_ADD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the product of two matrices plus a third, `A * B + C`, where at
 * least one argument has `var` scalars. The whole expression is one node on
 * the reverse-mode stack. The argument order follows `fma(x, y, z) = x * y +
 * z`.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 * @tparam Mat3 type of the matrix or expression added to the product
 * @param A first matrix
 * @param B second matrix
 * @param C matrix added to the product
 * @return `A * B + C`
 * @throw std::invalid_argument if the number of columns of `A` does not
 * match the number of rows of `B`, or `C` does not have the dimensions of
 * `A * B`
 */
template <typename Mat1, typename Mat2, typename Mat3,
          require_all_eigen_t<Mat1, Mat2, Mat3>* = nullptr,
          require_any_st_var<Mat1, Mat2, Mat3>* = nullptr,
          require_all_not_complex_t<value_type_t<Mat1>, value_type_t<Mat2>,
                                    value_type_t<Mat3>>* = nullptr>
inline auto multiply_add(Mat1&& A, Mat2&& B, Mat3&& C) {
  check_multiplicable("multiply_add", "A", A, "B", B);
  check_size_match("multiply_add", "Rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("multiply_add", "Columns of B", B.cols(), "columns of C",
                   C.cols());
  arena_t<Mat1> arena_A(std::forward<Mat1>(A));
  arena_t<Mat2> arena_B(std::forward<Mat2>(B));
  arena_t<Mat3> arena_C(std::forward<Mat3>(C));
  using return_t
      = return_var_matrix_t<decltype(arena_A * arena_B), Mat1, Mat2, Mat3>;
  arena_t<return_t> res
      = value_of(arena_A) * value_of(arena_B) + value_of(arena_C);
  reverse_pass_callback([arena_A, arena_B, arena_C, res]() mutable {
    // Eigen's matrix product cannot read the adjoint view directly
    const auto res_adj = res.adj().eval();
    if constexpr (is_autodiff_v<Mat1>) {
      arena_A.adj() += res_adj * value_of(arena_B).transpose();
    }
    if constexpr (is_autodiff_v<Mat2>) {
      arena_B.adj() += value_of(arena_A).transpose() * res_adj;
    }
    if constexpr (is_autodiff_v<Mat3>) {
      arena_C.adj() += res_adj;
    }
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
