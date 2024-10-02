#ifndef STAN_MATH_FWD_CONSTRAINT_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_FWD_CONSTRAINT_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/tcrossprod.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/constraint/unit_vector_constrain.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename EigMat,
          require_eigen_col_vector_vt<is_fvar, EigMat>* = nullptr>
inline auto unit_vector_constrain(const EigMat& y) {
  using eig_partial = partials_type_t<value_type_t<EigMat>>;
  promote_scalar_t<eig_partial, EigMat> y_val(value_of(y));
  plain_type_t<EigMat> unit_vector_y(y_val.size());
  unit_vector_y.val() = unit_vector_constrain(y_val);

  eig_partial squared_norm = dot_self(y_val);
  eig_partial norm = sqrt(squared_norm);
  eig_partial inv_norm = inv(norm);
  Eigen::Matrix<eig_partial, Eigen::Dynamic, Eigen::Dynamic> J
      = divide(tcrossprod(y_val), -norm * squared_norm);

  for (Eigen::Index m = 0; m < y_val.size(); ++m) {
    J.coeffRef(m, m) += inv_norm;
    for (Eigen::Index k = 0; k < y_val.size(); ++k) {
      unit_vector_y.coeffRef(k).d_ = J.coeff(k, m);
    }
  }
  return unit_vector_y;
}

template <typename EigMat, typename T,
          require_eigen_vt<is_fvar, EigMat>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
inline auto unit_vector_constrain(const EigMat& y, T& lp) {
  const auto& y_ref = to_ref(y);
  const value_type_t<EigMat> squared_norm = dot_self(y_ref);
  lp -= 0.5 * squared_norm;
  return unit_vector_constrain(y_ref);
}

}  // namespace math
}  // namespace stan
#endif
