#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>
#include <stan/math/rev/core.hpp>

#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaled step size, the function
 * tolerance, and the maximum number of steps.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] dat Continuous data vector for the equation system.
 * @param[in] dat_int Integer data vector for the equation system.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] scaling_step_size Scaled-step stopping tolerance. If
 *            a Newton step is smaller than the scaling step
 *            tolerance, the code breaks, assuming the solver is no
 *            longer making significant progress (i.e. is stuck)
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 *  * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T>
Eigen::VectorXd algebra_solver_newton(
    const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
    const Eigen::VectorXd& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    double scaling_step_size = 1e-3, double function_tolerance = 1e-6,
    long int max_num_steps = 200) {  // NOLINT(runtime/int)
  algebra_solver_check(x, y, dat, dat_int, function_tolerance, max_num_steps);
  check_nonnegative("algebra_solver", "scaling_step_size", scaling_step_size);

  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       value_of(f(x, y, dat, dat_int, msgs)),
                       "the vector of unknowns, x,", x);

  return kinsol_solve(f, value_of(x), y, dat, dat_int, 0, scaling_step_size,
                      function_tolerance, max_num_steps);
}

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaled step size, the function
 * tolerance, and the maximum number of steps.
 *
 * Overload the previous definition to handle the case where y
 * is a vector of parameters (var). The overload calls the
 * algebraic solver defined above and builds a vari object on
 * top, using the algebra_solver_vari class.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] dat Continuous data vector for the equation system.
 * @param[in] dat_int Integer data vector for the equation system.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] scaling_step_size Scaled-step stopping tolerance. If
 *            a Newton step is smaller than the scaling step
 *            tolerance, the code breaks, assuming the solver is no
 *            longer making significant progress (i.e. is stuck)
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @return theta Vector of solutions to the system of equations.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double scaling_step_size = 1e-3,
    double function_tolerance = 1e-6,
    long int max_num_steps = 200) {  // NOLINT(runtime/int)

  Eigen::VectorXd theta_dbl = algebra_solver_newton(
      f, x, value_of(y), dat, dat_int, msgs, scaling_step_size,
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
