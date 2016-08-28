#ifndef STAN_MATH_FWD_MAT_FUN_DOGLEG_HPP
#define STAN_MATH_FWD_MAT_FUN_DOGLEG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <stan/math/prim/mat/fun/inverse.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/mat/fun/to_fvar.hpp>

namespace stan {
  namespace math {

    template<typename T, typename F1, typename F2>
    inline
    Eigen::Matrix<fvar<T>, Eigen::Dynamic, 1>
    dogleg(const Eigen::Matrix<fvar<T>, Eigen::Dynamic, 1>& x,
           const F1, const F2) {

      Eigen::Vector<T, Eigen::Dynamic, 1> theta(x.rows());
      for (i = 0 i < x.rows(), i++) theta.coeffRef(i) = x.coeff(i).val_;
      Eigen::Vector<T, Eigen::Dynamic, 1> theta_star =
        stan::math::dogleg(theta, F1, F2);
      return to_fvar(theta_star, stan::math::inverse(F2(theta_star));
    }

  }
}
#endif
