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
 * @tparam EigMat1 type of the first matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A first matrix (m x n)
 * @param B second matrix (n x m)
 * @return trace of A * B
 * @throw std::invalid_argument if A and B have incompatible dimensions
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<std::is_arithmetic, EigMat1, EigMat2>* = nullptr>
inline return_type_t<EigMat1, EigMat2> trace_dot(const EigMat1& A,
                                                 const EigMat2& B) {
  check_size_match("trace_dot", "A.cols()", A.cols(), "B.rows()", B.rows());
  check_size_match("trace_dot", "A.rows()", A.rows(), "B.cols()", B.cols());
  return A.cwiseProduct(B.transpose()).sum();
}

}  // namespace math
}  // namespace stan

#endif
