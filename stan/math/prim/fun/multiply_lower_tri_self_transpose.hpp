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
 * @throw std::domain_error If the input matrix is not square.
 */
template <typename EigMat, typename = require_eigen_t<EigMat>>
inline auto multiply_lower_tri_self_transpose(EigMat&& L) {
  int K = L.rows();
  if (K == 0) {
    return Eigen::Matrix<value_type_t<EigMat>, -1, -1>(0, 0);
  }
  if (K == 1) {
    Eigen::Matrix<value_type_t<EigMat>, -1, -1> result(1, 1);
    result(0) = square(L(0));  // first elt, so don't need double idx
    return result;
  }
  return (L.template triangularView<Eigen::Lower>() * L.transpose()).eval();
}

}  // namespace math
}  // namespace stan

#endif
