#ifndef STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/to_matrix.hpp>
#include <Eigen/dense>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Return analytical form, which will be used to compute Jacobians.
 */
template <typename F1, typename F2, typename F3, typename F4, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
quadratic_optimizer_analytical
  (const F1& fh,
   const F2& fv,
   const F3& fa,
   const F4& fb,
   const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
   const Eigen::VectorXd& delta,
   const Eigen::VectorXd& x) {

  typedef T scalar;

  using Eigen::Matrix;
  using Eigen::Dynamic;

  double tol = 1e-10;  // FIX ME - add this as an argument.

  Matrix<scalar, Dynamic, Dynamic> H = fh(theta, delta);
  Matrix<scalar, Dynamic, Dynamic> v = fv(theta, delta);
  Matrix<scalar, Dynamic, 1> A = fa(theta, delta);
  scalar b = fb(theta, delta);
  int size = x.size();

  Matrix<T, Dynamic, 1> x_analytical(x.size());

  Matrix<scalar, Dynamic, 1> sum_x_H(size);
  for (int j = 0; j < size; j++) {
    Matrix<scalar, Dynamic, 1> H_row = H.row(j);
    sum_x_H(j) = dot_product(H_row, x) - x(j) * H(j, j);
  }

  scalar sum_A_H = 0;
  for (int j = 0; j < size; j++) if (x(j) > tol)
    sum_A_H += A(j) * A(j) / H(j, j);

  scalar sum_lambda = 0;
  for (int k = 0; k < size; k++) {
    if (x(k) > tol)
      sum_lambda += A(k) / H(k, k) * (v(k) + 0.5 * sum_x_H(k));
  }
  // Here we write the equality constraint A.x + b = 0,
  // as opposed to A.x = b, hence the substraction.
  sum_lambda -= b;

  Matrix<scalar, Dynamic, 1> x_an(size);
  for (int i = 0; i < x.size(); i++)
    x_an(i) = (1 / H(i, i)) * (A(i) * sum_lambda / sum_A_H
                                 - v(i) - 0.5 * sum_x_H(i));
  return x_an;
}

/**
* A functor which can be passed to the jacobian function.
*/
template <typename Fh, typename Fv, typename Fa, typename Fb>
struct f_theta {
  Eigen::VectorXd delta_;
  Eigen::VectorXd x_;
  Fh fh_;
  Fv fv_;
  Fa fa_;
  Fb fb_;

  f_theta(const Fh& fh,
          const Fv& fv,
          const Fa& fa,
          const Fb& fb,
          const Eigen::VectorXd& delta,
          const Eigen::VectorXd& x) : delta_(delta), x_(x),
          fh_(fh), fv_(fv), fa_(fa), fb_(fb) {}

  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  inline operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    return quadratic_optimizer_analytical(fh_, fv_, fa_, fb_,
                                          theta, delta_, x_);
  }
};

/**
 * The vari class for the quadratic optimizer.
 */
template <typename F, typename T>
struct quadratic_optimizer_vari : public vari {
  /** vector of parameters */
  vari** theta_;
  /** number of parameters */
  int theta_size_;
  /** vector of solution */
  vari** x_;
  /** number of unknowns */
  int x_size_;
  /** functor for jacobian function */
  F f_;
  /** Jacobian of the solution w.r.t the parameters */
  double* Jx_theta_;

  quadratic_optimizer_vari(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
                           const Eigen::VectorXd& x,
                           const F& f)
    : vari(x(0)),
      theta_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
          theta.size())),
      theta_size_(theta.size()),
      x_(ChainableStack::instance().memalloc_.alloc_array<vari*>(x.size())),
      x_size_(x.size()),
      f_(f),
      Jx_theta_(ChainableStack::instance().memalloc_.alloc_array<double>(
        x_size_ * theta_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;

    for (int i = 0; i < theta_size_; i++) theta_[i] = theta(i).vi_;

    x_[0] = this;
    for (int i = 1; i < x_size_; i++) x_[i] = new vari(x(i), false);

    // Compute the Jacobian matrix and store in array.
    Eigen::VectorXd f_eval;  // dummy argument required by jacobian function.
    MatrixXd J;
    jacobian(f_, value_of(theta), f_eval, J);

    Map<MatrixXd>(&Jx_theta_[0], x_size_, theta_size_) = J;
  }

  void chain() {
    for (int j = 0; j < theta_size_; j++)
      for (int i = 0; i < x_size_; i++)
        theta_[i]->adj_ += x_[i]->adj_ * Jx_theta_[j * x_size_ + i];
  }
};

/**
 * Return the solution to the specified quadratic optimization
 * problem. Evaluation function with theta a parameter of doubles.
 */
template <typename Fh, typename Fv, typename Fa, typename Fb>
Eigen::VectorXd
quadratic_optimizer(const Fh& fh, 
                    const Fv& fv, 
                    const Fa& fa, 
                    const Fb& fb,
                    const Eigen::VectorXd& theta,
                    const Eigen::VectorXd& delta,
                    int n) {
  using Eigen::VectorXd;

  VectorXd x;
  double f_value;  // declare here to remove warning message

  VectorXd b(1);
  b(0) = fb(theta, delta);

  VectorXd A_v = fa(theta, delta);
  Eigen::MatrixXd A(1, A_v.size());
  A = to_matrix(A_v);

  f_value = Eigen::solve_quadprog(fh(theta, delta),
                                  fv(theta, delta), A, b,
                                  Eigen::MatrixXd::Identity(n, n),
                                  Eigen::VectorXd::Zero(n),
                                  x);
  return x;
}

/**
 * Return the solution to the specified quadratic optimization
 * problem. Evaluation of the function with theta a parameter of var.
 */
template <typename Fh, typename Fv, typename Fa, typename Fb, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
quadratic_optimizer(const Fh& fh,
                    const Fv& fv,
                    const Fa& fa,
                    const Fb& fb,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
                    const Eigen::VectorXd& delta,
                    int n) {
  Eigen::VectorXd
    x_dbl = quadratic_optimizer(fh, fv, fa, fb, value_of(theta), delta, n);

  // Construct vari
  typedef f_theta<Fh, Fv, Fa, Fb> f_vari;
  quadratic_optimizer_vari<f_vari, T>* vi0
    = new quadratic_optimizer_vari<f_vari, T>(theta, x_dbl,
      f_vari(fh, fv, fa, fb, delta, x_dbl));

  Eigen::Matrix<var, Eigen::Dynamic, 1> x(x_dbl.size());
  x(0) = var(vi0->x_[0]);
  for (int i = 1; i < x.size(); i++)
    x(i) = var(vi0->x_[i]);

  return x;
}

}  // namespace math
}  // namespace stan

#endif
