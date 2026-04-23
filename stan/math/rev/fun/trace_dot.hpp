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
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param A first matrix (m x n)
 * @param B second matrix (n x m)
 * @return trace of A * B
 * @throw std::invalid_argument if A and B have incompatible dimensions
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
inline var trace_dot(const Mat1& A, const Mat2& B) {
  check_size_match("trace_dot", "A.cols()", A.cols(), "B.rows()", B.rows());
  check_size_match("trace_dot", "A.rows()", A.rows(), "B.cols()", B.cols());

  var res;

  if constexpr (is_autodiff_v<Mat1> && is_autodiff_v<Mat2>) {
    arena_t<promote_scalar_t<var, Mat1>> arena_A = A;
    arena_t<promote_scalar_t<var, Mat2>> arena_B = B;

    res = value_of(arena_A).cwiseProduct(value_of(arena_B).transpose()).sum();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += res.adj() * value_of(arena_B).transpose();
      } else {
        arena_A.adj() += res.adj() * value_of(arena_B).transpose();
      }
      if constexpr (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += res.adj() * value_of(arena_A).transpose();
      } else {
        arena_B.adj() += res.adj() * value_of(arena_A).transpose();
      }
    });
  } else if constexpr (is_autodiff_v<Mat2>) {
    arena_t<promote_scalar_t<double, Mat1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, Mat2>> arena_B = B;

    res = arena_A.cwiseProduct(value_of(arena_B).transpose()).sum();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if constexpr (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += res.adj() * arena_A.transpose();
      } else {
        arena_B.adj() += res.adj() * arena_A.transpose();
      }
    });
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_A = A;
    arena_t<promote_scalar_t<double, Mat2>> arena_B = value_of(B);

    res = value_of(arena_A).cwiseProduct(arena_B.transpose()).sum();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += res.adj() * arena_B.transpose();
      } else {
        arena_A.adj() += res.adj() * arena_B.transpose();
      }
    });
  }

  return res;
}

}  // namespace math
}  // namespace stan
#endif
