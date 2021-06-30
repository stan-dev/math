#ifndef STAN_MATH_PRIM_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP
#define STAN_MATH_PRIM_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Private interface for calling the Powell solver. Users should call the Powell
 * solver through `algebra_solver_powell` or `algebra_solver_powell_impl`.
 *
 * @tparam F type of equation system function, curried with respect to inputs
 * @tparam T type of elements in the x vector
 * @tparam Args types of additional parameters to the equation system functor
 *
 * @param[in] f Functor that evaluates the system of equations, curried with
 *            respect to its input (i.e., "y") values.
 * @param[in] x Vector of starting values (initial guess).
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance determines the convergence criteria
 *            for the solution.
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @return theta_dbl Double vector of solutions to the system of equations.
 * @pre x has size greater than zero.
 * @pre x has only finite elements.
 * @pre f returns finite values when passed any value of x and the given args.
 * @pre relative_tolerance is non-negative.
 * @pre function_tolerance is non-negative.
 * @pre max_num_steps is positive.
 * @throw <code>std::domain_error</code> solver exceeds max_num_steps.
 * @throw <code>std::domain_error</code> if the norm of the solution exceeds
 * the function tolerance.
 */
template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T>* = nullptr>
T& algebra_solver_powell_call_solver(
    const F& f, T& x, std::ostream* const msgs,
    const double relative_tolerance, const double function_tolerance,
    const int64_t max_num_steps, const Args&... args) {
  // Construct the solver
  hybrj_functor_solver<decltype(f)> hfs(f);
  Eigen::HybridNonLinearSolver<decltype(hfs)> solver(hfs);

  // Compute theta_dbl
  solver.parameters.xtol = relative_tolerance;
  solver.parameters.maxfev = max_num_steps;
  solver.solve(x);

  // Check if the max number of steps has been exceeded
  if (solver.nfev >= max_num_steps) {
    [&]() STAN_COLD_PATH {
    throw_domain_error("algebra_solver", "maximum number of iterations",
                       max_num_steps, "(", ") was exceeded in the solve.");
    }();
  }

  // Check solution is a root
  double system_norm = f(x).stableNorm();
  if (system_norm > function_tolerance) {
    [&]() STAN_COLD_PATH {
    std::ostringstream message;
    message << "the norm of the algebraic function is " << system_norm
            << " but should be lower than the function "
            << "tolerance:";
    throw_domain_error("algebra_solver", message.str().c_str(),
                       function_tolerance, "",
                       ". Consider decreasing the relative tolerance and "
                       "increasing max_num_steps.");
    }();
  }

  return x;
}

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * Use Powell's dogleg solver.
 *
 * The user can also specify the relative tolerance
 * (xtol in Eigen's code), the function tolerance,
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * This function is overloaded to handle both constant and var-type parameters.
 * This overload handles non-var parameters, and checks the input and calls the
 * algebraic solver only.
 *
 * @tparam F type of equation system function
 * @tparam T type of elements in the x vector
 * @tparam Args types of additional parameters to the equation system functor
 *
 * @param[in] f Functor that evaluates the system of equations.
 * @param[in] x Vector of starting values (initial guess).
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance determines the convergence criteria
 *            for the solution.
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps maximum number of function evaluations.
 * @param[in] args additional parameters to the equation system functor.
 * @return theta Vector of solutions to the system of equations.
 * @pre f returns finite values when passed any value of x and the given args.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * elements.
 * @throw <code>std::invalid_argument</code> if relative_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>std::domain_error</code> solver exceeds max_num_steps.
 * @throw <code>std::domain_error</code> if the norm of the solution exceeds
 * the function tolerance.
 */
template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_all_st_arithmetic<Args...>* = nullptr>
Eigen::VectorXd algebra_solver_powell_impl(const F& f, const T& x,
                                           std::ostream* const msgs,
                                           const double relative_tolerance,
                                           const double function_tolerance,
                                           const int64_t max_num_steps,
                                           const Args&... args) {
  const auto& x_ref = to_ref(x);
  auto x_val = to_ref(value_of(x_ref));

  // Curry the input function w.r.t. y
  auto f_wrt_x = [&f, msgs, &args...](const auto& x) { return f(x, msgs, args...); };

  check_nonzero_size("algebra_solver_powell", "initial guess", x_val);
  check_finite("algebra_solver_powell", "initial guess", x_val);
  check_nonnegative("alegbra_solver_powell", "relative_tolerance",
                    relative_tolerance);
  check_nonnegative("algebra_solver_powell", "function_tolerance",
                    function_tolerance);
  check_positive("algebra_solver_powell", "max_num_steps", max_num_steps);
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       f_wrt_x(x_ref), "the vector of unknowns, x,", x_ref);

  // Solve the system
  return algebra_solver_powell_call_solver(f_wrt_x, x_val, msgs,
                                            relative_tolerance,
                                            function_tolerance, max_num_steps);
}

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * Use Powell's dogleg solver.
 *
 * The user can also specify the relative tolerance
 * (xtol in Eigen's code), the function tolerance,
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * Signature to maintain backward compatibility, will be removed
 * in the future.
 *
 * @tparam F type of equation system function
 * @tparam T1 type of elements in the x vector
 * @tparam T2 type of elements in the y vector
 *
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
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite
 * elements.
 * @throw <code>std::invalid_argument</code> if relative_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>std::domain_error</code> solver exceeds max_num_steps.
 * @throw <code>std::domain_error</code> if the norm of the solution exceeds
 * the function tolerance.
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2>* = nullptr>
Eigen::Matrix<value_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_powell(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* const msgs = nullptr,
    const double relative_tolerance = 1e-10,
    const double function_tolerance = 1e-6,
    const int64_t max_num_steps = 1e+3) {
  return algebra_solver_powell_impl(algebra_solver_adapter<F>(f), x, msgs,
                                    relative_tolerance, function_tolerance,
                                    max_num_steps, y, dat, dat_int);
}

/**
 * Return the solution to the specified system of algebraic
 * equations given an initial guess, and parameters and data,
 * which get passed into the algebraic system.
 * Use Powell's dogleg solver.
 *
 * The user can also specify the relative tolerance
 * (xtol in Eigen's code), the function tolerance,
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * Signature to maintain backward compatibility, will be removed
 * in the future.
 *
 * @tparam F type of equation system function
 * @tparam T1 type of elements in the x vector
 * @tparam T2 type of elements in the y vector
 *
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
 * @pre f returns finite values when passed any value of x and the given y, dat,
 *        and dat_int.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite
 * elements.
 * @throw <code>std::invalid_argument</code> if relative_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::domain_error</code>) if solver exceeds max_num_steps.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::domain_error</code>) if the norm of the solution exceeds the
 * function tolerance.
 */
template <typename F, typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2>* = nullptr>
Eigen::Matrix<value_type_t<T2>, Eigen::Dynamic, 1> algebra_solver(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    const double relative_tolerance = 1e-10,
    const double function_tolerance = 1e-6,
    const int64_t max_num_steps = 1e+3) {
  return algebra_solver_powell(f, x, y, dat, dat_int, msgs, relative_tolerance,
                               function_tolerance, max_num_steps);
}

}  // namespace math
}  // namespace stan

#endif
