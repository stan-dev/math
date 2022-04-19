#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/kinsol_solve.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
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
 * This overload handles non-autodiff parameters.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @tparam Args types of additional parameters to the equation system functor
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] scaling_step_size Scaled-step stopping tolerance. If
 *            a Newton step is smaller than the scaling step
 *            tolerance, the code breaks, assuming the solver is no
 *            longer making significant progress (i.e. is stuck)
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @param[in] args Additional parameters to the equation system functor.
 * @return theta Vector of solutions to the system of equations.
 * @pre f returns finite values when passed any value of x and the given args.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>std::domain_error if solver exceeds max_num_steps.
 */
template <typename F, typename T, typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_all_st_arithmetic<Args...>* = nullptr>
Eigen::VectorXd algebra_solver_newton_impl(const F& f, const T& x,
                                           std::ostream* const msgs,
                                           const double scaling_step_size,
                                           const double function_tolerance,
                                           const int64_t max_num_steps,
                                           const Args&... args) {
  const auto& x_ref = to_ref(value_of(x));

  check_nonzero_size("algebra_solver_newton", "initial guess", x_ref);
  check_finite("algebra_solver_newton", "initial guess", x_ref);
  check_nonnegative("algebra_solver_newton", "scaling_step_size",
                    scaling_step_size);
  check_nonnegative("algebra_solver_newton", "function_tolerance",
                    function_tolerance);
  check_positive("algebra_solver_newton", "max_num_steps", max_num_steps);

  return kinsol_solve(f, x_ref, scaling_step_size, function_tolerance,
                      max_num_steps, 1, 10, KIN_LINESEARCH, msgs, args...);
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
 * This overload handles var parameters.
 *
 * The Jacobian \(J_{xy}\) (i.e., Jacobian of unknown \(x\) w.r.t. the parameter
 * \(y\)) is calculated given the solution as follows. Since
 * \[
 *   f(x, y) = 0,
 * \]
 * we have (\(J_{pq}\) being the Jacobian matrix \(\tfrac {dq} {dq}\))
 * \[
 *   - J_{fx} J_{xy} = J_{fy},
 * \]
 * and therefore \(J_{xy}\) can be solved from system
 * \[
 *  - J_{fx} J_{xy} = J_{fy}.
 * \]
 * Let \(eta\) be the adjoint with respect to \(x\); then to calculate
 * \[
 *   \eta J_{xy},
 * \]
 * we solve
 * \[
 *   - (\eta J_{fx}^{-1}) J_{fy}.
 * \]
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @tparam Args types of additional parameters to the equation system functor
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] scaling_step_size Scaled-step stopping tolerance. If
 *            a Newton step is smaller than the scaling step
 *            tolerance, the code breaks, assuming the solver is no
 *            longer making significant progress (i.e. is stuck)
 * @param[in] function_tolerance determines whether roots are acceptable.
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @param[in] args Additional parameters to the equation system functor.
 * @return theta Vector of solutions to the system of equations.
 * @pre f returns finite values when passed any value of x and the given args.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>std::domain_error if solver exceeds max_num_steps.
 */
template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T>* = nullptr,
          require_any_st_var<T_Args...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> algebra_solver_newton_impl(
    const F& f, const T& x, std::ostream* const msgs,
    const double scaling_step_size, const double function_tolerance,
    const int64_t max_num_steps, const T_Args&... args) {
  const auto& x_ref = to_ref(value_of(x));
  auto arena_args_tuple = make_chainable_ptr(std::make_tuple(eval(args)...));
  auto args_vals_tuple = apply(
      [&](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      *arena_args_tuple);

  check_nonzero_size("algebra_solver_newton", "initial guess", x_ref);
  check_finite("algebra_solver_newton", "initial guess", x_ref);
  check_nonnegative("algebra_solver_newton", "scaling_step_size",
                    scaling_step_size);
  check_nonnegative("algebra_solver_newton", "function_tolerance",
                    function_tolerance);
  check_positive("algebra_solver_newton", "max_num_steps", max_num_steps);

  // Solve the system
  Eigen::VectorXd theta_dbl = apply(
      [&](const auto&... vals) {
        return kinsol_solve(f, x_ref, scaling_step_size, function_tolerance,
                            max_num_steps, 1, 10, KIN_LINESEARCH, msgs,
                            vals...);
      },
      args_vals_tuple);

  auto f_wrt_x = [&](const auto& x) {
    return apply([&](const auto&... args) { return f(x, msgs, args...); },
                 args_vals_tuple);
  };

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(f_wrt_x, theta_dbl, f_x, Jf_x);

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  arena_t<ret_type> ret = theta_dbl;
  auto Jf_x_T_lu_ptr
      = make_unsafe_chainable_ptr(Jf_x.transpose().partialPivLu());  // Lu

  reverse_pass_callback(
      [f, ret, arena_args_tuple, Jf_x_T_lu_ptr, msgs]() mutable {
        Eigen::VectorXd eta = -Jf_x_T_lu_ptr->solve(ret.adj().eval());

        // Contract with Jacobian of f with respect to y using a nested reverse
        // autodiff pass.
        {
          nested_rev_autodiff rev;

          Eigen::VectorXd ret_val = ret.val();
          auto x_nrad_ = apply(
              [&ret_val, &f, msgs](const auto&... args) {
                return eval(f(ret_val, msgs, args...));
              },
              *arena_args_tuple);
          x_nrad_.adj() = eta;
          grad();
        }
      });

  return ret_type(ret);
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
 * Signature to maintain backward compatibility, will be removed
 * in the future.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system.
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
          require_all_eigen_vector_t<T1, T2>* = nullptr>
Eigen::Matrix<scalar_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* const msgs = nullptr,
    const double scaling_step_size = 1e-3,
    const double function_tolerance = 1e-6,
    const long int max_num_steps = 200) {  // NOLINT(runtime/int)
  return algebra_solver_newton_impl(algebra_solver_adapter<F>(f), x, msgs,
                                    scaling_step_size, function_tolerance,
                                    max_num_steps, y, dat, dat_int);
}

}  // namespace math
}  // namespace stan

#endif
