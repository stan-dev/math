#ifndef STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
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
  Eigen::VectorXd x;
  double f_value;  // declare here to remove warning message
  
  f_value = Eigen::solve_quadprog(fh(theta, delta),
                                  fv(theta, delta),
                                  fa(theta, delta),
                                  fb(theta, delta),
                                  Eigen::MatrixXd::Identity(n, n),
                                  Eigen::VectorXd::Zero(n),
                                  x);
  return x;
}

}  // namespace math
}  // namespace stan

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

   typedef typename T scalar;
  
   using Eigen::Matrix;
   using Eigen::Dynamic;
   using stan::math::sum;
   
   double tol = 1e-10;  // FIX ME - add this as an argument.

   Matrix<scalar, Dynamic, Dynamic> H = fh(theta, delta);
   Matrix<scalar, Dynamic, Dynamic> v = fv(theta, delta);
   Matrix<scalar, Dynamic, 1> A = fa(theta, delta);
   scalar b = fb(theta, delta);
   int size = x.size();

   Matrix<T, Dynamic, 1> x_analytical(x.size());

   Matrix<scalar, Dynamic, 1> sum_x_H;
   for (int j = 0; j < size; j++)
     sum_x_H(j) = dot_product(H.row(i), x) - x(i) * H(i, i);

   scalar sum_A_H = 0;
   for (int j = 0; j < size; j++) if (x(j) > tol) sum_A_H += A(j)^2 / H(j, j);

   scalar sum_lambda = 0;
   for (int k = 0; k < size; k++) {
     if (x(k) > tol) 
       sum_lambda += A(k) / H(k, k) * (v(k) + 0.5 * sum_x_H(k));
   }
   sum_lambda += b;

   for (int i = 0; i < x.size(); i++)
     x(i) = (1 / H(i, i)) * (A(i) * sum_lambda / sum_A_H
       + v(i) - 0.5 * sum_x_H(i));
}






#endif
