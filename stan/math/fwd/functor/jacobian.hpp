#ifndef STAN_MATH_FWD_FUNCTOR_JACOBIAN_HPP
#define STAN_MATH_FWD_FUNCTOR_JACOBIAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T, typename F>
void jacobian(const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
              Eigen::Matrix<T, Eigen::Dynamic, 1>& fx,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& J) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<fvar<T>, Dynamic, 1> x_fvar(x.size());
  J.resize(x_fvar.size(), x.size());
  fx.resize(x_fvar.size());
  for (int k = 0; k < x.size(); ++k) {
    x_fvar(k) = fvar<T>(x(k), 0);
  }
  x_fvar(0) = fvar<T>(x(0), 1);
  Matrix<fvar<T>, Dynamic, 1> fx_fvar = f(x_fvar);
  fx = fx_fvar.val();
  J.col(0) = fx_fvar.d();
  const fvar<T> switch_fvar(0, 1);  // flips the tangents on and off
  for (int i = 1; i < x.size(); ++i) {
    x_fvar(i - 1) -= switch_fvar;
    x_fvar(i) += switch_fvar;
    Matrix<fvar<T>, Dynamic, 1> fx_fvar = f(x_fvar);
    J.col(i) = fx_fvar.d();
  }
}

}  // namespace math
}  // namespace stan
#endif
