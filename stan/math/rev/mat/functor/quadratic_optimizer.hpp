#ifndef STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
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
 * Return the solution to the specified quadratic optimization
 * problem. Evaluation function with theta a parameter of doubles.
 */
template <typename F1, typename F2, typename F3, typename F4>
Eigen::VectorXd
quadratic_optimizer(const F1& fh, 
                    const F2& fv, 
                    const F3& fa, 
                    const F4& fb,
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

}  // namespace math
}  // namespace stan

#endif
