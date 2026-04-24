#ifndef STAN_MATH_PRIM_FUN_TRACE_DOT_HPP
#define STAN_MATH_PRIM_FUN_TRACE_DOT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Compute the trace of the product of two matrices,
 * \f$ \text{tr}(A \cdot B) = \sum_{i,j} A_{ij} B_{ji} \f$.
 *
 * This is more efficient than computing the full product and
 * taking the trace, as it avoids forming the intermediate matrix.
 *
 * @tparam EigMat1 A type either inheriting from `Eigen::DenseBase` or a
 * `var_value` with an inner type inheriting from `Eigen::DenseBase`
 * @tparam EigMat2 A type either inheriting from `Eigen::DenseBase` or a
 * `var_value` with an inner type inheriting from `Eigen::DenseBase`
 *
 * @param A first matrix (m x n)
 * @param B second matrix (n x m)
 * @return trace of A * B
 * @throw std::invalid_argument if A and B have incompatible dimensions
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<std::is_arithmetic, EigMat1, EigMat2>* = nullptr>
inline auto trace_dot(EigMat1&& A, EigMat2&& B) {
  check_size_match("trace_dot", "A.cols()", A.cols(), "B.rows()", B.rows());
  check_size_match("trace_dot", "A.rows()", A.rows(), "B.cols()", B.cols());
  return make_holder([](auto&& A_, auto&& B_) {
    return A_.cwiseProduct(B_.transpose()).sum();
  }, std::forward<EigMat1>(A), std::forward<EigMat2>(B));
}

}  // namespace math
}  // namespace stan

#endif
