#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>

#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * 
 * Use Kinsol's Newton solver.
 */
template <typename F, typename T>
Eigen::VectorXd algebra_solver_newton(
  const F&f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
  const Eigen::VectorXd& y, const std::vector<double>& dat,
  const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
  double relative_tolerance = 1e-10, double function_tolerance = 1e-6,
  long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)

  return kinsol_solve(f, value_of(x), y, dat, dat_int, 0,
                      function_tolerance, max_num_steps, 1e-3);
  }

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * 
 * Overload for case where theta is a vector of parameters.
 * 
 * Use Kinsol's Newton solver.
 */
template <typename F, typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-10,
    double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)

  Eigen::VectorXd theta_dbl
    = algebra_solver_newton(f, x, value_of(y), dat, dat_int,
                            msgs, relative_tolerance,
                            function_tolerance, max_num_steps);

  typedef system_functor<F, double, double, false> Fy;
  typedef system_functor<F, double, double, true> Fs;
  typedef hybrj_functor_solver<Fs, F, double, double> Fx;
  Fx fx(Fs(), f, value_of(x), value_of(y), dat, dat_int, msgs);

  // Construct vari
  algebra_solver_vari<Fy, F, T2, Fx>* vi0
    = new algebra_solver_vari<Fy, F, T2, Fx>(Fy(), f, value_of(x), y, dat,
                                             dat_int, theta_dbl, fx, msgs);
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(x.size());
  theta(0) = var(vi0->theta_[0]);
  for (int i = 1; i < x.size(); ++i)
    theta(i) = var(vi0->theta_[i]);

  return theta;
}

}  // namespace math
}  // namespace stan

#endif
