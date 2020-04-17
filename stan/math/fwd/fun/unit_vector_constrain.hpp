#ifndef STAN_MATH_FWD_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_FWD_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/tcrossprod.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/unit_vector_constrain.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline auto unit_vector_constrain(EigMat&& y) {
  using Eigen::Matrix;
  using std::sqrt;
  using eig_mat = std::decay_t<EigMat>;
  using ref_inner = typename eig_mat::PlainObject;
  using eig_index = index_type_t<EigMat>;
  using eig_partial = partials_type_t<value_type_t<EigMat>>;
  using eig_value = value_type_t<EigMat>;
  using partial_mat = Eigen::Matrix<eig_partial, eig_mat::RowsAtCompileTime,
                                    eig_mat::ColsAtCompileTime>;
  partial_mat y_t(y.rows(), y.cols());
  const Eigen::Ref<const ref_inner, Eigen::Aligned16, Eigen::Stride<0, 0>>&
      y_ref
      = y;
  for (eig_index k = 0; k < y.size(); ++k) {
    y_t.coeffRef(k) = y_ref.coeff(k).val_;
  }

  partial_mat unit_vector_y_t = unit_vector_constrain(y_t);
  ref_inner unit_vector_y(y.size());
  for (eig_index k = 0; k < y.size(); ++k) {
    unit_vector_y.coeffRef(k).val_ = unit_vector_y_t.coeff(k);
  }

  eig_partial squared_norm = dot_self(y_t);
  eig_partial norm = sqrt(squared_norm);
  eig_partial inv_norm = inv(norm);
  Eigen::Matrix<eig_partial, Eigen::Dynamic, Eigen::Dynamic> J
      = divide(tcrossprod(y_t), -norm * squared_norm);

  for (eig_index m = 0; m < y.size(); ++m) {
    J.coeffRef(m, m) += inv_norm;
    for (eig_index k = 0; k < y.size(); ++k) {
      unit_vector_y.coeffRef(k).d_ = J.coeff(k, m);
    }
  }
  return unit_vector_y;
}

template <typename EigMat, require_eigen_vt<is_fvar, EigMat>* = nullptr>
inline auto unit_vector_constrain(EigMat&& y, value_type_t<EigMat>& lp) {
  const value_type_t<EigMat> squared_norm = dot_self(y);
  lp -= 0.5 * squared_norm;
  return unit_vector_constrain(std::forward<EigMat>(y));
}

}  // namespace math
}  // namespace stan
#endif
