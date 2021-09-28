#ifndef STAN_MATH_FWD_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_FWD_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::RowsAtCompileTime>
multiply_lower_tri_self_transpose(const EigMat& m) {
  if (m.rows() == 0) {
    return {};
  }
  plain_type_t<EigMat> L = m.template triangularView<Eigen::Lower>();

  return multiply(L, L.transpose());
}

}  // namespace math
}  // namespace stan
#endif
