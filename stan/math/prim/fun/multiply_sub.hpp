#ifndef STAN_MATH_PRIM_FUN_MULTIPLY_SUB_HPP
#define STAN_MATH_PRIM_FUN_MULTIPLY_SUB_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a matrix minus the product of two matrices, `C - A * B`. The
 * argument order follows `multiply_add(A, B, C) = A * B + C`.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 * @tparam Mat3 type of the matrix or expression the product is subtracted
 * from
 * @param A first matrix
 * @param B second matrix
 * @param C matrix the product is subtracted from
 * @return `C - A * B`
 * @throw std::invalid_argument if the number of columns of `A` does not
 * match the number of rows of `B`, or `C` does not have the dimensions of
 * `A * B`
 */
template <typename Mat1, typename Mat2, typename Mat3,
          require_all_eigen_t<Mat1, Mat2, Mat3>* = nullptr,
          require_all_st_arithmetic<Mat1, Mat2, Mat3>* = nullptr>
inline auto multiply_sub(Mat1&& A, Mat2&& B, Mat3&& C) {
  check_multiplicable("multiply_sub", "A", A, "B", B);
  check_size_match("multiply_sub", "Rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("multiply_sub", "Columns of B", B.cols(), "columns of C",
                   C.cols());
  return make_holder(
      [](auto&& A_, auto&& B_, auto&& C_) {
        return std::forward<decltype(C_)>(C_)
               - std::forward<decltype(A_)>(A_)
                     * std::forward<decltype(B_)>(B_);
      },
      std::forward<Mat1>(A), std::forward<Mat2>(B), std::forward<Mat3>(C));
}

}  // namespace math
}  // namespace stan

#endif
