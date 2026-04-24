#ifndef STAN_MATH_FWD_FUN_TRACE_DOT_HPP
#define STAN_MATH_FWD_FUN_TRACE_DOT_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/trace.hpp>

namespace stan {
namespace math {

/**
 * Compute the trace of the product of two matrices with
 * forward-mode autodiff support.
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
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_vt_fvar<EigMat1, EigMat2>* = nullptr>
inline return_type_t<EigMat1, EigMat2> trace_dot(EigMat1&& A,
                                                 EigMat2&& B) {
  check_size_match("trace_dot", "A.cols()", A.cols(), "B.rows()", B.rows());
  check_size_match("trace_dot", "A.rows()", A.rows(), "B.cols()", B.cols());
  return trace(multiply(std::forward<EigMat1>(A), std::forward<EigMat2>(B)));
}

}  // namespace math
}  // namespace stan
#endif
