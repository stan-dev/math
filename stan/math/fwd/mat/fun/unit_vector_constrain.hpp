#ifndef STAN_MATH_FWD_MAT_FUN_UNIT_VECTOR_CONSTRAIN_HPP
#define STAN_MATH_FWD_MAT_FUN_UNIT_VECTOR_CONSTRAIN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/mat/fun/dot_self.hpp>
#include <stan/math/fwd/scal/fun/sqrt.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/tcrossprod.hpp>
#include <stan/math/prim/mat/fun/unit_vector_constrain.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>

namespace stan {
  namespace math {

    template <typename T, int R, int C>
    inline Eigen::Matrix<fvar<T>, R, C>
    unit_vector_constrain(const Eigen::Matrix<fvar<T>, R, C>& y) {
      using Eigen::Matrix;
      using stan::math::tcrossprod;
      using stan::math::dot_self;
      using stan::math::unit_vector_constrain;

      Matrix<T, R, C> y_t(y.size());
      for (int k = 0; k < y.size(); ++k)
        y_t.coeffRef(k) = y.coeff(k).val_;

      Matrix<T, R, C> unit_vector_y_t
        = unit_vector_constrain(y_t);
      Matrix<fvar<T>, R, C> unit_vector_y(y.size());
      for (int k = 0; k < y.size(); ++k)
        unit_vector_y.coeffRef(k).val_ = unit_vector_y_t.coeff(k);

      const T squared_norm = dot_self(y_t);
      const T norm = sqrt(squared_norm);
      const T inv_norm = inv(norm);
      Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J
        = tcrossprod(y_t) / (-norm * squared_norm);

      // for each input position
      for (int m = 0; m < y.size(); ++m) {
        J.coeffRef(m, m) += inv_norm;
        // for each output position
        for (int k = 0; k < y.size(); ++k) {
          // chain from input to output
          unit_vector_y.coeffRef(k).d_ = J.coeff(k, m);
        }
      }
      return unit_vector_y;
    }

    template <typename T, int R, int C>
    inline Eigen::Matrix<fvar<T>, R, C>
    unit_vector_constrain(const Eigen::Matrix<fvar<T>, R, C>& y, fvar<T>& lp) {
      using stan::math::dot_self;
      const fvar<T> squared_norm = dot_self(y);
      lp -= 0.5 * squared_norm;
      return unit_vector_constrain(y);
    }

  }
}
#endif
