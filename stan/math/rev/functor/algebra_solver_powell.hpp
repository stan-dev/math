#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
#include <stan/math/prim/functor/algebra_solver_powell.hpp>
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
 * Use Powell's dogleg solver.
 *
 * The user can also specify the relative tolerance
 * (xtol in Eigen's code), the function tolerance,
 * and the maximum number of steps (maxfev in Eigen's code).
 *
 * This function is overloaded to handle both constant and var-type parameters.
 * This overload handles var parameters, and checks the input, calls the
 * algebraic solver, and appropriately handles derivative propagation through
 * the `reverse_pass_callback`.
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
 * @param[in] max_num_steps  maximum number of function evaluations.
 * @param[in] args Additional parameters to the equation system functor.
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
template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T>* = nullptr,
          require_any_st_var<T_Args...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> algebra_solver_powell_impl(
    const F& f, const T& x, std::ostream* const msgs,
    const double relative_tolerance, const double function_tolerance,
    const int64_t max_num_steps, const T_Args&... args) {
  const auto& x_ref = to_ref(x);
  auto x_val = to_ref(value_of(x_ref));
  auto arena_args_tuple = std::make_tuple(to_arena(args)...);
  auto args_vals_tuple = apply(
      [&](const auto&... args) {
        return std::make_tuple(to_ref(value_of(args))...);
      },
      arena_args_tuple);

  // Curry the input function w.r.t. y
  auto f_wrt_x = [&args_vals_tuple, &f, msgs](const auto& x) {
    return apply([&x, &f, msgs](const auto&... args) { return f(x, msgs, args...); },
                 args_vals_tuple);
  };

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
  algebra_solver_powell_call_solver(
      f_wrt_x, x_val, msgs, relative_tolerance, function_tolerance,
      max_num_steps);

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(f_wrt_x, x_val, f_x, Jf_x);

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  auto Jf_xT_lu_ptr
      = make_unsafe_chainable_ptr(Jf_x.transpose().partialPivLu());  // Lu

  arena_t<ret_type> ret = x_val;

  reverse_pass_callback([f, ret, arena_args_tuple, Jf_xT_lu_ptr,
                         msgs]() mutable {
    // Contract specificities with inverse Jacobian of f with respect to x.
    Eigen::VectorXd eta = -Jf_xT_lu_ptr->solve(ret.adj().eval());

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      nested_rev_autodiff rev;
      Eigen::VectorXd ret_val = ret.val();
      auto x_nrad_ = apply(
          [&](const auto&... args) { return eval(f(ret_val, msgs, args...)); },
          arena_args_tuple);
      x_nrad_.adj() = eta;
      grad();
    }
  });

  return ret_type(ret);
}

}  // namespace math
}  // namespace stan

#endif
