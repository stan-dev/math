#ifndef STAN_MATH_REV_FUNCTOR_JACOBIAN_HPP
#define STAN_MATH_REV_FUNCTOR_JACOBIAN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

template <typename F>
void jacobian(const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
              Eigen::Matrix<double, Eigen::Dynamic, 1>& fx,
              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& J) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  // Run nested autodiff in this scope
  nested_rev_autodiff nested;

  Matrix<var, Dynamic, 1> x_var(x);
  Matrix<var, Dynamic, 1> fx_var = f(x_var);
  fx.resize(fx_var.size());
  J.resize(x.size(), fx_var.size());
  fx = fx_var.val();
  grad(fx_var(0).vi_);
  J.col(0) = x_var.adj();
  for (int i = 1; i < fx_var.size(); ++i) {
    nested.set_zero_all_adjoints();
    grad(fx_var(i).vi_);
    J.col(i) = x_var.adj();
  }
  J.transposeInPlace();
}

}  // namespace math
}  // namespace stan
#endif
