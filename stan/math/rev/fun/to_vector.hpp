#ifndef STAN_MATH_REV_FUN_TO_VECTOR_HPP
#define STAN_MATH_REV_FUN_TO_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Reshape a `var_value<Matrix>` to a `var_value<ColumnVector>`.
 * @tparam EigMat Inner type of the `var_value` that must inherit from
 *  `Eigen::EigenBase`.
 * @param x A var whose inner matrix type is to be reshaped to an Column matrix.
 * @return A view of the original `x` with inner value and adjoint matrices
 *  mapped to a column vector.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline auto to_vector(const var_value<EigMat>& x) {
  using view_type = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>;
  return var_value<view_type>(new vari_view<view_type>(
      view_type(x.vi_->val_.data(), x.rows() * x.cols()),
      view_type(x.vi_->adj_.data(), x.rows() * x.cols())));
}

}  // namespace math
}  // namespace stan
#endif
