#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>

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
 * Use a Newton solver.
 *
 * The user can also specify the relative tolerance 
 * (xtol in Eigen's code), the function tolerance, 
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * Throw an exception if the norm of f(x), where f is the
 * output of the algebraic system and x the proposed solution,
 * is greater than the function tolerance. We here use the
 * norm as a metric to measure how far we are from the origin (0).
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @param[in] f Functor that evaluates the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] dat continuous data vector for the equation system.
 * @param[in] dat_int integer data vector for the equation system.
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance determines the convergence criteria
 *            for the solution.
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @return theta Vector of solutions to the system of equations.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if relative_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if the norm of the solution exceeds the
 * function tolerance.
 */
template <typename F, typename T>
Eigen::VectorXd algebra_solver_newton(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
    const Eigen::VectorXd& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    double relative_tolerance = 1e-10, double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  algebra_solver_check(x, y, dat, dat_int, 
                       relative_tolerance, function_tolerance, max_num_steps);

  Eigen::VectorXd x_dbl = value_of(x);
  Eigen::VectorXd fx;
  Eigen::MatrixXd J;
  Eigen::MatrixXd direction;

  system_functor<F, T, double, true> system(f, x, y, dat, dat_int, msgs);

  for (int i = 0; i <= max_num_steps; i++) {
    if (i == max_num_steps) {
      std::ostringstream message;
      message << "algebra_solver_newton: max number of iterations: "
              << max_num_steps << " exceeded.";
      throw boost::math::evaluation_error(message.str());
    }

    jacobian(system, x_dbl, fx, J);
    direction = - mdivide_left(J, fx);
    x_dbl += direction;  // CHECK - add linesearch method?

    // Check if the solution is acceptable
    if (fx.norm() <= function_tolerance) break;
  }

  // Check dimension unknowns equals dimension of system output
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       fx, "the vector of unknowns, x,",
                       x);

  // CHECK - do we need additional checks? Does the function break if
  // we don't reach the max number of steps?

  return x_dbl;
}

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * Use a Newton solver.
 *
 * The user can also specify the relative tolerance 
 * (xtol in Eigen's code), the function tolerance, 
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * Overload the previous definition to handle the case where y
 * is a vector of parameters (var). The overload calls the
 * algebraic solver defined above and builds a vari object on
 * top, using the algebra_solver_vari class.
 *
 * @tparam F type of equation system function.
 * @tparam T1  Type of elements in x vector.
 * @tparam T2  Type of elements in y vector.
 * @param[in] f Functor that evaluates the system of equations.
 * @param[in] x Vector of starting values (initial guess).
 * @param[in] y parameter vector for the equation system.
 * @param[in] dat continuous data vector for the equation system.
 * @param[in] dat_int integer data vector for the equation system.
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance determines the convergence criteria
 *            for the solution.
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @return theta Vector of solutions to the system of equations.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if relative_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if the norm of the solution exceeds the
 * function tolerance.
 */
template <typename F, typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-10,
    double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)

  bool test_kinsol = 1;

  Eigen::VectorXd theta_dbl;

  if (!test_kinsol) {
    theta_dbl
      = algebra_solver_newton(f, x, value_of(y), dat, dat_int, 0,
                              relative_tolerance,
                              function_tolerance, max_num_steps);
  } else {
    theta_dbl
      = kinsol_solve(f, x, value_of(y), dat, dat_int, 0,
                     function_tolerance, max_num_steps);
  }

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
