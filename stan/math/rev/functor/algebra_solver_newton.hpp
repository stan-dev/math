#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_NEWTON_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/functor/kinsol_solve.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/** Implementation of ordinary newton solver. */
template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T>* = nullptr,
          require_all_st_arithmetic<T_Args...>* = nullptr>
Eigen::VectorXd algebra_solver_newton_impl(
    const F& f, const T& x, std::ostream* msgs, double scaling_step_size,
    double function_tolerance, long int max_num_steps, const Eigen::VectorXd& y,
    const T_Args&... args) {  // NOLINT(runtime/int)
  const auto& x_val = to_ref(value_of(x));
  auto args_vals_tuple = std::make_tuple(y, to_ref(args)...);

  check_nonzero_size("algebra_solver_newton", "initial guess", x_val);
  check_finite("algebra_solver_newton", "initial guess", x_val);
  check_nonnegative("algebra_solver_newton", "scaling_step_size", scaling_step_size);
  check_nonnegative("algebra_solver_newton", "function_tolerance", function_tolerance);
  check_positive("algebra_solver_newton", "max_num_steps", max_num_steps);

  return kinsol_solve(f, x_val, scaling_step_size,
                      function_tolerance, max_num_steps, 1, 10, KIN_LINESEARCH,
                      msgs, y, args...);
}

/** Implementation of autodiff newton solver. */
template <typename F, typename T, typename... T_Args,
          require_eigen_vector_t<T>* = nullptr,
          require_any_st_var<T_Args...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> algebra_solver_newton_impl(
    const F& f, const T& x, std::ostream* msgs, double scaling_step_size,
    double function_tolerance, long int max_num_steps,
    const T_Args&... args) {  // NOLINT(runtime/int)
  const auto& x_val = to_ref(value_of(x));
  auto arena_args_tuple = std::make_tuple(to_arena(args)...);
  auto args_vals_tuple = std::make_tuple(to_ref(value_of(args))...);

  check_nonzero_size("algebra_solver_newton", "initial guess", x_val);
  check_finite("algebra_solver_newton", "initial guess", x_val);
  check_nonnegative("algebra_solver_newton", "scaling_step_size", scaling_step_size);
  check_nonnegative("algebra_solver_newton", "function_tolerance", function_tolerance);
  check_positive("algebra_solver_newton", "max_num_steps", max_num_steps);

  // Solve the system
  Eigen::VectorXd theta_dbl = apply(
      [&](const auto&... vals) {
        return kinsol_solve(f, x_val, scaling_step_size,
                            function_tolerance, max_num_steps, 1, 10,
                            KIN_LINESEARCH, msgs, vals...);
      },
      args_vals_tuple);

  // Evaluate and store the Jacobian.
  auto myfunc = [&](const auto& x) {
    return apply([&](const auto&... args) { return f(x, msgs, args...); },
                 args_vals_tuple);
  };

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(myfunc, theta_dbl, f_x, Jf_x);

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  auto arena_Jf_x = to_arena(Jf_x);

  arena_t<ret_type> ret = theta_dbl;

  reverse_pass_callback([f, ret, arena_args_tuple, arena_Jf_x, msgs]() mutable {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Contract specificities with inverse Jacobian of f with respect to x.
    VectorXd ret_adj = ret.adj();
    VectorXd eta = -arena_Jf_x.transpose().fullPivLu().solve(ret_adj);

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      nested_rev_autodiff rev;

      VectorXd ret_val = ret.val();
      auto x_nrad_ = apply(
          [&](const auto&... args) { return eval(f(ret_val, msgs, args...)); },
          arena_args_tuple);
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
          require_all_eigen_vector_t<T1, T2>* = nullptr>
Eigen::Matrix<scalar_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_newton(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    double scaling_step_size = 1e-3, double function_tolerance = 1e-6,
    long int max_num_steps = 200) {  // NOLINT(runtime/int)
  return algebra_solver_newton_impl(algebra_solver_adapter<F>(f), x, msgs,
                                    scaling_step_size, function_tolerance,
                                    max_num_steps, y, dat, dat_int);
}

}  // namespace math
}  // namespace stan

#endif
