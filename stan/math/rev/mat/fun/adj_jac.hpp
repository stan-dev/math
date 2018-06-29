#ifndef STAN_MATH_REV_MAT_FUN_ADJ_JAC_HPP
#define STAN_MATH_REV_MAT_FUN_ADJ_JAC_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace stan {
namespace math {

/*
 * The adjoint Jacobian vari is a special vari for wrapping up user defined ops that supply an apply operation (that performs an operation on some input) and a multiply_adjoint_jacobian call (that computes the result of the product of the transpose of the adjoint times the Jacobian of the operator).
 *
 * The adjoint Jacobian vari itself does nothing other than connect the supplied operation to the autodiff stack and coordinate movement of data between all the involved varis so that someone can implement autodiff functions without ever having to deal with autodiff types.
 */
template <typename F>
struct adj_jac_vari : public vari {
  F* f_;
  int N_;
  vari** x_vi_;
  int M_;
  vari** y_vi_;

  adj_jac_vari(F* f, const Eigen::Matrix<var, Eigen::Dynamic, 1>& x, int M, vari** y_vi)
    : vari(0.0), f_(f),
        N_(x.size()), x_vi_(ChainableStack::instance().memalloc_.alloc_array<vari*>(N_)),
	M_(M), y_vi_(y_vi) {
    for (int n = 0; n < N_; ++n)
      x_vi_[n] = x(n).vi_;
  }

  void chain() {
    Eigen::Matrix<double, Eigen::Dynamic, 1> y_adj(M_);
    for (int m = 0; m < M_; ++m)
      y_adj(m) = y_vi_[m] -> adj_;
    Eigen::Matrix<double, Eigen::Dynamic, 1> y_adj_jac = f_ -> multiply_adjoint_jacobian(y_adj);
    for (int n = 0; n < N_; ++n)
      x_vi_[n] -> adj_ += y_adj_jac(n);
  }
};

/*
 * adj_jac_apply makes it possible to write reasonably efficient reverse-mode autodiff code without ever touching Stan's autodiff internals
 *
 * Mathematically, to use a function in reverse mode autodiff, you need to be able to evaluate the function (y = f(x)) and the Jacobian of that function (df(x)/dx).
 *
 * Pretend there exists some large, complicated function, L(x), which contains our smaller function f(x). The goal of autodiff is to compute dL/dx. If we break the large function into pieces:
 *
 * y = f(x)
 * L = g(y)
 *
 * Then if we were given dL/dy, then we could compute dL/dx by the product dL/dy * dy/dx
 *
 * Because y = f(x), dy/dx is just df(x)/dx, the Jacobian of the function we're trying to define. In vector form,
 *
 * dL/dx = (dL/dy)' * df(x)/dx
 *
 * So implementing f(x) and the product above is all that is required mathematically to implement reverse-mode autodiff for a function.
 *
 * adj_jac_apply takes as a template argument a functor that supplies:
 *   Eigen::VectorXd apply(const Eigen::VectorXd& x) and,
 *   Eigen::VectorXd multiply_adjoint_jacobian(const Eigen::VectorXd& adj)
 *
 *  apply is f(x) and multiply_adjoint_jacobian is responsible for computing the adjoint transpose Jacobian product (which frequently does not require the calculation of the full Jacobian).
 *
 * The functor supplied to adj_jac_apply must be careful to allocate itself as well as any other variables it defines on the autodiff memory stack because its destructor will never be called! (and so memory will leak if allocated anywhere else)
 *
 * @tparam F functor to be connected to the autodiff stack
 * @param x input to the functor
 * @return the result of the specified operation wrapped up in vars
 */
template <class F>
Eigen::Matrix<var, Eigen::Dynamic, 1> adj_jac_apply(const Eigen::Matrix<var, Eigen::Dynamic, 1>& x) {
  F *f = new F();
  Eigen::Matrix<double, Eigen::Dynamic, 1> val_y = f -> apply(value_of(x));

  int M = val_y.size();
  vari** y_vi = ChainableStack::instance().memalloc_.alloc_array<vari*>(M);
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(M);
  for (int m = 0; m < M; ++m) {
    y_vi[m] = new vari(val_y(m), false);
    y(m) = var(y_vi[m]);
  }

  new adj_jac_vari<F>(f, x, M, y_vi);
  return y;
}

}  // namespace math
}  // namespace stan
#endif
