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
  using eigen_ref
      = Eigen::Ref<const ref_inner, Eigen::Aligned16, Eigen::Stride<0, 0>>;
  using eig_index = index_type_t<EigMat>;
  using eig_partial = partials_type_t<value_type_t<EigMat>>;
  using eig_value = value_type_t<EigMat>;
  using partial_mat = Eigen::Matrix<eig_partial, eig_mat::RowsAtCompileTime,
                                    eig_mat::ColsAtCompileTime>;
  partial_mat y_val(y.rows(), y.cols());
  const eigen_ref& y_ref = y;
  y_val = value_of(y_ref);

  partial_mat unit_vector_y_val = unit_vector_constrain(y_val);
  ref_inner unit_vector_y(y.size());
  unit_vector_y.val() = unit_vector_y_val;

  eig_partial squared_norm = dot_self(y_val);
  eig_partial norm = sqrt(squared_norm);
  eig_partial inv_norm = inv(norm);
  Eigen::Matrix<eig_partial, Eigen::Dynamic, Eigen::Dynamic> J
      = divide(tcrossprod(y_val), -norm * squared_norm);

  for (eig_index m = 0; m < y.size(); ++m) {
    J.coeffRef(m, m) += inv_norm;
    for (eig_index k = 0; k < y.size(); ++k) {
      unit_vector_y.coeffRef(k).d_ = J.coeff(k, m);
    }
  }
  return unit_vector_y;
}

template <typename EigMat, typename T,
          require_eigen_vt<is_fvar, EigMat>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
inline auto unit_vector_constrain(EigMat&& y, T& lp) {
  const value_type_t<EigMat> squared_norm = dot_self(y);
  lp -= 0.5 * squared_norm;
  return unit_vector_constrain(std::forward<EigMat>(y));
}

}  // namespace math
}  // namespace stan
#endif
