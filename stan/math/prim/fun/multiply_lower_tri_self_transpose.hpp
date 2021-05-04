#ifndef STAN_MATH_PRIM_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_PRIM_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Returns the result of multiplying the lower triangular
 * portion of the input matrix by its own transpose.
 *
 * @param L Matrix to multiply.
 * @return The lower triangular values in L times their own
 * transpose.
 */
template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
          require_not_st_autodiff<EigMat>* = nullptr>
inline matrix_d multiply_lower_tri_self_transpose(const EigMat& L) {
  int K = L.rows();
  if (K == 0) {
    return L;
  }
  if (K == 1) {
    matrix_d result(1, 1);
    result.coeffRef(0) = square(L.coeff(0, 0));
    return result;
  }
  int J = L.cols();
  matrix_d LLt(K, K);
  matrix_d Lt = L.transpose();
  for (int m = 0; m < K; ++m) {
    int k = (J < m + 1) ? J : m + 1;
    LLt(m, m) = Lt.col(m).head(k).squaredNorm();
    for (int n = (m + 1); n < K; ++n) {
      LLt(n, m) = LLt(m, n) = Lt.col(m).head(k).dot(Lt.col(n).head(k));
    }
  }
  return LLt;
}

}  // namespace math
}  // namespace stan

#endif
