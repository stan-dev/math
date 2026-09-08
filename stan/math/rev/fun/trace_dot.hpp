#ifndef STAN_MATH_REV_FUN_TRACE_DOT_HPP
#define STAN_MATH_REV_FUN_TRACE_DOT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/trace_dot.hpp>

namespace stan {
namespace math {

/**
 * Compute the trace of the product of two matrices with autodiff support.
 *
 * \f$ \text{tr}(A \cdot B) = \sum_{i,j} A_{ij} B_{ji} \f$
 *
 * The gradients are:
 * \f$ \frac{\partial}{\partial A} \text{tr}(A B) = B^T \f$,
 * \f$ \frac{\partial}{\partial B} \text{tr}(A B) = A^T \f$.
 *
 * @tparam Mat1 A type either inheriting from `Eigen::DenseBase` or a
 * `var_value` with an inner type inheriting from `Eigen::DenseBase`
 * @tparam Mat2 A type either inheriting from `Eigen::DenseBase` or a
 * `var_value` with an inner type inheriting from `Eigen::DenseBase`
 *
 * @param A first matrix (m x n)
 * @param B second matrix (n x m)
 * @return trace of A * B
 * @throw std::invalid_argument if A and B have incompatible dimensions
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
inline var trace_dot(Mat1&& A, Mat2&& B) {
  check_size_match("trace_dot", "A.cols()", A.cols(), "B.rows()", B.rows());
  check_size_match("trace_dot", "A.rows()", A.rows(), "B.cols()", B.cols());
  if constexpr (is_autodiff_v<Mat1> && is_autodiff_v<Mat2>) {
    arena_t<Mat1> arena_A(std::forward<Mat1>(A));
    arena_t<Mat2> arena_B(std::forward<Mat2>(B));
    auto res_val = arena_A.val().cwiseProduct(arena_B.val().transpose()).sum();
    return make_callback_var(res_val, [arena_A, arena_B](auto&& res) mutable {
      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += res.adj() * arena_B.val().transpose();
      } else {
        arena_A.adj() += res.adj() * arena_B.val().transpose();
      }
      if constexpr (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += res.adj() * arena_A.val().transpose();
      } else {
        arena_B.adj() += res.adj() * arena_A.val().transpose();
      }
    });
  } else if constexpr (is_autodiff_v<Mat2>) {
    arena_t<Mat1> arena_A(std::forward<Mat1>(A));
    arena_t<Mat2> arena_B(std::forward<Mat2>(B));
    auto res_val = arena_A.cwiseProduct(arena_B.val().transpose()).sum();
    return make_callback_var(res_val, [arena_A, arena_B](auto&& res) mutable {
      if constexpr (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += res.adj() * arena_A.transpose();
      } else {
        arena_B.adj() += res.adj() * arena_A.transpose();
      }
    });
  } else {
    arena_t<Mat1> arena_A(std::forward<Mat1>(A));
    arena_t<Mat2> arena_B(std::forward<Mat2>(B));
    auto res_val = arena_A.val().cwiseProduct(arena_B.transpose()).sum();
    return make_callback_var(res_val, [arena_A, arena_B](auto&& res) mutable {
      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += res.adj() * arena_B.transpose();
      } else {
        arena_A.adj() += res.adj() * arena_B.transpose();
      }
    });
  }
}

}  // namespace math
}  // namespace stan
#endif
