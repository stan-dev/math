#ifndef STAN_MATH_MIX_FUNCTOR_GRAD_TR_MAT_TIMES_HESSIAN_HPP
#define STAN_MATH_MIX_FUNCTOR_GRAD_TR_MAT_TIMES_HESSIAN_HPP

#include <stan/math/mix/functor/gradient_dot_vector.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

template <typename F>
void grad_tr_mat_times_hessian(
    const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_tr_MH) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  // Run nested autodiff in this scope
  nested_rev_autodiff nested;

  grad_tr_MH.resize(x.size());

  Matrix<var, Dynamic, 1> x_var(x.size());
  for (int i = 0; i < x.size(); ++i) {
    x_var(i) = x(i);
  }

  Matrix<fvar<var>, Dynamic, 1> x_fvar(x.size());

  var sum(0.0);
  Matrix<double, Dynamic, 1> M_n(x.size());
  for (int n = 0; n < x.size(); ++n) {
    for (int k = 0; k < x.size(); ++k) {
      M_n(k) = M(n, k);
    }
    for (int k = 0; k < x.size(); ++k) {
      x_fvar(k) = fvar<var>(x_var(k), k == n);
    }
    fvar<var> fx;
    fvar<var> grad_fx_dot_v;
    gradient_dot_vector<fvar<var>, double>(f, x_fvar, M_n, fx, grad_fx_dot_v);
    sum += grad_fx_dot_v.d_;
  }

  grad(sum.vi_);
  for (int i = 0; i < x.size(); ++i) {
    grad_tr_MH(i) = x_var(i).adj();
  }
}

}  // namespace math
}  // namespace stan
#endif
