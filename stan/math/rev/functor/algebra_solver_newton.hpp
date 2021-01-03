#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/functor/kinsol_solve.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
 *
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
 * @throw <code>std::domain_error</code> if solver exceeds max_num_steps.
 */
template <typename F, typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::VectorXd algebra_solver_newton(
    const F& f, const T& x, const Eigen::VectorXd& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double scaling_step_size = 1e-3,
    double function_tolerance = 1e-6,
    long int max_num_steps = 200) {  // NOLINT(runtime/int)
  const auto& x_eval = x.eval();
  algebra_solver_check(x_eval, y, dat, dat_int, function_tolerance,
                       max_num_steps);
  check_nonnegative("algebra_solver", "scaling_step_size", scaling_step_size);

  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       value_of(f(x_eval, y, dat, dat_int, msgs)),
                       "the vector of unknowns, x,", x);

  return kinsol_solve(f, value_of(x_eval), y, dat, dat_int, 0,
                      scaling_step_size, function_tolerance, max_num_steps);
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
 *
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
 * @throw <code>std::domain_error if solver exceeds max_num_steps.
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2>* = nullptr,
          require_st_var<T2>* = nullptr>
Eigen::Matrix<scalar_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    double scaling_step_size = 1e-3, double function_tolerance = 1e-6,
    long int max_num_steps = 200) {  // NOLINT(runtime/int)

  const auto& x_eval = x.eval();
  const auto& y_eval = y.eval();
  Eigen::VectorXd theta_dbl = algebra_solver_newton(
      f, x_eval, value_of(y_eval), dat, dat_int, msgs, scaling_step_size,
      function_tolerance, max_num_steps);

  typedef system_functor<F, double, double, false> Fy;
  typedef system_functor<F, double, double, true> Fs;
  typedef hybrj_functor_solver<Fs, F, double, double> Fx;
  Fx fx(Fs(), f, value_of(x_eval), value_of(y_eval), dat, dat_int, msgs);

  // Construct vari
  auto* vi0 = new algebra_solver_vari<Fy, F, scalar_type_t<T2>, Fx>(
      Fy(), f, value_of(x_eval), y_eval, dat, dat_int, theta_dbl, fx, msgs);
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(x.size());
  theta(0) = var(vi0->theta_[0]);
  for (int i = 1; i < x.size(); ++i)
    theta(i) = var(vi0->theta_[i]);

  return theta;
}

}  // namespace math
}  // namespace stan

#endif
