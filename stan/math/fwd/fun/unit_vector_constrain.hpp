#ifndef STAN_MATH_FWD_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_FWD_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/divide.hpp>
#include <stan/math/fwd/fun/dot_self.hpp>
#include <stan/math/fwd/fun/tcrossprod.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/prim/fun/divide.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/unit_vector_constrain.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> unit_vector_constrain(
    const Eigen::Matrix<fvar<T>, R, C>& y) {
  using Eigen::Matrix;
  using std::sqrt;

  Matrix<T, R, C> y_t(y.size());
  for (int k = 0; k < y.size(); ++k) {
    y_t.coeffRef(k) = y.coeff(k).val_;
  }

  Matrix<T, R, C> unit_vector_y_t = unit_vector_constrain(y_t);
  Matrix<fvar<T>, R, C> unit_vector_y(y.size());
  for (int k = 0; k < y.size(); ++k) {
    unit_vector_y.coeffRef(k).val_ = unit_vector_y_t.coeff(k);
  }

  T squared_norm = dot_self(y_t);
  T norm = sqrt(squared_norm);
  T inv_norm = inv(norm);
  Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J
      = divide(tcrossprod(y_t), -norm * squared_norm);

  for (int m = 0; m < y.size(); ++m) {
    J.coeffRef(m, m) += inv_norm;
    for (int k = 0; k < y.size(); ++k) {
      unit_vector_y.coeffRef(k).d_ = J.coeff(k, m);
    }
  }
  return unit_vector_y;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> unit_vector_constrain(
    const Eigen::Matrix<fvar<T>, R, C>& y, fvar<T>& lp) {
  fvar<T> squared_norm = dot_self(y);
  lp -= 0.5 * squared_norm;
  return unit_vector_constrain(y);
}

}  // namespace math
}  // namespace stan
#endif
