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
  
  // Eigen::VectorXd test_constrain(2);
  // test_constrain << -1, -10;
  
  f_value = Eigen::solve_quadprog(fh(theta, delta),
                                  fv(theta, delta),
                                  fa(theta, delta),
                                  fb(theta, delta),
                                  Eigen::MatrixXd::Identity(n, n),
                                  // test_constrain,
                                  Eigen::VectorXd::Zero(n),
                                  x);
  // std::cout << "f_value: " << f_value << "\n";
  return x;
}

}  // namespace math
}  // namespace stan

#endif
