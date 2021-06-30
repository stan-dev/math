#ifndef STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP
#define STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/kinsol_data.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
#include <stan/math/prim/functor/algebra_solver_fp.hpp>
#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a fixed pointer to the specified system of algebraic
 * equations of form
 * \[
 *   x = F(x; theta)
 * \]
 * given an initial guess \(x\), and parameters \(theta\) and data. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaling controls, the function
 * tolerance, and the maximum number of steps.
 *
 * This function is overloaded to handle both constant and var-type parameters.
 * This overload handles var parameters, and checks the input, calls the
 * algebraic solver, and appropriately handles derivative propagation through
 * the `reverse_pass_callback`.
 *
 * The Jacobian \(J_{xy}\) (i.e., Jacobian of unknown \(x\) w.r.t. the parameter
 * \(y\)) is calculated given the solution as follows. Since
 * \[
 *   x - f(x, y) = 0,
 * \]
 * we have (\(J_{pq}\) being the Jacobian matrix \(\tfrac {dq} {dq}\))
 * \[
 *   J_{xy} - J_{fx} J_{xy} = J_{fy},
 * \]
 * and therefore \(J_{xy}\) can be solved from system
 * \[
 *   (I - J_{fx}) * J_{xy} = J_{fy}.
 * \]
 * Let \(eta\) be the adjoint with respect to \(x\); then to calculate
 * \[
 *   \eta J_{xy},
 * \]
 * we solve
 * \[
 *   (\eta * (I - J_{fx})^{-1}) * J_{fy}.
 * \]
 * (This is virtually identical to the Powell and Newton solvers, except
 * \(-J_{fx}\) has been replaced by \((I - J_{fx}\).)
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 *
 * @param[in] f functor that evaluated the system of equations.
 * @param[in] x vector of starting values.
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] u_scale diagonal scaling matrix elements \(Du\)
 *                    such that \(Du x\) has all components roughly the same
 *                    magnitude when \(x\) is close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_scale diagonal scaling matrix elements such
 *                    that \(Df (x - f(x))\) has all components roughly the same
 *                    magnitude when \(x\) is not too close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] function_tolerance Function-norm stopping tolerance.
 * @param[in] max_num_steps maximum number of function evaluations.
 * @param[in] args additional parameters to the equation system functor.
 * @pre f returns finite values when passed any value of x and the given args.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T, typename T_u, typename T_f, typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_any_st_var<Args...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> algebra_solver_fp_impl(
    const F& f, const T& x, std::ostream* const msgs,
    const std::vector<T_u>& u_scale, const std::vector<T_f>& f_scale,
    const double function_tolerance, const int max_num_steps,
    const Args&... args) {
  auto arena_args_tuple = std::make_tuple(to_arena(args)...);

  auto args_vals_tuple = apply(
      [](const auto&... args) { return std::make_tuple(value_of(args)...); },
      arena_args_tuple);

  auto f_wrt_x = [&msgs, &f, &args_vals_tuple](const auto& x) {
    return apply(
        [&x, &msgs, &f](const auto&... args) { return f(x, msgs, args...); },
        args_vals_tuple);
  };

  // FP solution
  Eigen::VectorXd theta_dbl = apply(
      [&f, function_tolerance, &u_scale, &f_scale, &msgs,
        max_num_steps, x_val = value_of(x)](const auto&... vals) {
        return kinsol_solve_fp(f, x_val, function_tolerance,
                               max_num_steps, u_scale, f_scale, msgs, vals...);
      },
      args_vals_tuple);

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(f_wrt_x, theta_dbl, f_x, Jf_x);

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  arena_t<ret_type> ret = theta_dbl;

  auto Jf_xT_lu_ptr
      = make_unsafe_chainable_ptr((Eigen::MatrixXd::Identity(x.size(), x.size()) - Jf_x).transpose().eval().partialPivLu());  // Lu

  reverse_pass_callback([f, ret, arena_args_tuple, Jf_xT_lu_ptr, msgs]() mutable {
    // Contract specificities with inverse Jacobian of f with respect to x.
    Eigen::VectorXd eta = Jf_xT_lu_ptr->solve(ret.adj().eval());

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      nested_rev_autodiff rev;
      auto x_nrad_ = apply(
          [&](const auto&... args) { return eval(f(ret.val(), msgs, args...)); },
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
