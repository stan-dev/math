#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * The vari class for the algebraic solver. We compute the  Jacobian of
 * the solutions with respect to the parameters using the implicit
 * function theorem. The call to Jacobian() occurs outside the call to
 * chain() -- this prevents malloc issues.
 */
template <typename Fy, typename T, typename Fx>
struct algebra_solver_vari : public vari {
  /** number of parameters */
  int y_size_;
  /** value of parameters */
  arena_t<Eigen::Matrix<T, Eigen::Dynamic, 1>> y_val_;
  /** System functor (f_) w.r.t. inputs (y_) */
  Fy fy_;
  /** array of parameters */
  vari** y_;
  /** number of unknowns */
  int x_size_;
  /** value of unknowns */
  arena_t<Eigen::VectorXd&> x_val_;
  /** array of unknowns */
  vari** x_;
  /** Jacobian of f w.r.t. outputs (x_) */
  const Eigen::MatrixXd Jf_x_;

  algebra_solver_vari(Fy& fy, const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                      Fx& fx, const Eigen::VectorXd& x)
      : vari(x(0)),
        y_size_(y.size()),
        y_val_(y),
        fy_(fy),
        y_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(y_size_)),
        x_size_(x.size()),
        x_val_(x),
        Jf_x_(fx.get_jacobian(x_val_)),
        x_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(x_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;

    for (int i = 0; i < y_size_; ++i) {
      y_[i] = y(i).vi_;
    }

    x_[0] = this;
    for (int i = 1; i < x_size_; ++i) {
      x_[i] = new vari(x(i), false);
    }
  }

  void chain() {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Compute (transpose of) specificities with respect to x.
    VectorXd x_bar_(x_size_);
    for (int i = 0; i < x_size_; ++i) {
      x_bar_[i] = x_[i]->adj_;
    }

    // Contract specificities with inverse Jacobian of f with respect to x.
    VectorXd eta_ = -Jf_x_.transpose().fullPivLu().solve(x_bar_);

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      stan::math::nested_rev_autodiff rev;
      Matrix<var, Eigen::Dynamic, 1> y_nrad_ = y_val_;
      auto x_nrad_ = stan::math::eval(fy_(y_nrad_));
      x_nrad_.adj() = eta_;
      stan::math::grad();
    }
  }
};  // namespace math

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
 * Throw an exception if the norm of f(x), where f is the
 * output of the algebraic system and x the proposed solution,
 * is greater than the function tolerance. We here use the
 * norm as a metric to measure how far we are from the origin (0).
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 *
 * @param[in] f Functor that evaluates the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or as a
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
 * @throw <code>std::domain_error</code> solver exceeds max_num_steps.
 * @throw <code>std::domain_error</code> if the norm of the solution exceeds
 * the function tolerance.
 */
template <typename F, typename T, require_eigen_vector_t<T>* = nullptr>
Eigen::VectorXd algebra_solver_powell(
    const F& f, const T& x, const Eigen::VectorXd& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-10,
    double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  /*const auto& x_eval = x.eval();
  const auto& x_val = (value_of(x_eval)).eval();
  algebra_solver_check(x_val, y, dat, dat_int, function_tolerance,
                       max_num_steps);
  check_nonnegative("alegbra_solver", "relative_tolerance", relative_tolerance);

  // Create functor for algebraic system
  using Fs = system_functor<F, double, double, true>;
  using Fx = hybrj_functor_solver<Fs, F, double, double>;
  Fx fx(Fs(), f, x_val, y, dat, dat_int, msgs);
  Eigen::HybridNonLinearSolver<Fx> solver(fx);

  // Check dimension unknowns equals dimension of system output
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       fx.get_value(x_val), "the vector of unknowns, x,", x);

  // Solve the system
  return algebra_solver_powell_(solver, fx, x, y, dat, dat_int, msgs,
                                relative_tolerance, function_tolerance,
                                max_num_steps);*/
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
 * Overload the previous definition to handle the case where y
 * is a vector of parameters (var). The overload calls the
 * algebraic solver defined above and builds a vari object on
 * top, using the algebra_solver_vari class.
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
          require_all_eigen_vector_t<T1, T2>* = nullptr,
          require_st_var<T2>* = nullptr>
Eigen::Matrix<value_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_powell(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, std::ostream* msgs = nullptr,
    double relative_tolerance = 1e-10, double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  const auto& x_eval = x.eval();
  auto arena_y = to_arena(y);
  auto arena_dat = to_arena(dat);
  auto arena_dat_int = to_arena(dat_int);
  const auto& x_val = (value_of(x_eval)).eval();
  const auto& y_val = (value_of(arena_y)).eval();

  algebra_solver_check(x_val, y, dat, dat_int, function_tolerance,
                       max_num_steps);
  check_nonnegative("alegbra_solver", "relative_tolerance", relative_tolerance);

  // Construct the Powell solver

  auto myfunc = [&](const auto& x) {
    return f(x, y_val, dat, dat_int, msgs);
  };

  hybrj_functor_solver<decltype(myfunc)> fx(myfunc);
  Eigen::HybridNonLinearSolver<decltype(fx)> solver(fx);

  // Check dimension unknowns equals dimension of system output
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       fx.get_value(x_val), "the vector of unknowns, x,", x);

  // Solve the system
  Eigen::VectorXd theta_dbl = algebra_solver_powell_(
      solver, fx, x_eval, y_val, dat, dat_int, 0, relative_tolerance,
      function_tolerance, max_num_steps);

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(myfunc, theta_dbl, f_x, Jf_x);

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  auto arena_Jf_x = to_arena(Jf_x);

  arena_t<ret_type> ret = theta_dbl;

  reverse_pass_callback([f, ret, arena_y, arena_dat, arena_dat_int, arena_Jf_x, msgs]() mutable {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Contract specificities with inverse Jacobian of f with respect to x.
    std::cout << "ret_adj: " << ret.adj().transpose() << std::endl;
    std::cout << "Jfx: " << arena_Jf_x << std::endl;
    VectorXd ret_adj = ret.adj();
    VectorXd eta = -arena_Jf_x.transpose().fullPivLu().solve(ret_adj);

    std::cout << "eta: " << eta.transpose() << std::endl;

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      stan::math::nested_rev_autodiff rev;
      Matrix<var, Eigen::Dynamic, 1> y_nrad_ = arena_y.val();

      std::cout << "y_val: " << arena_y.val().transpose() << std::endl;
      VectorXd ret_val = ret.val();
      std::cout << "ret_val: " << ret_val.transpose() << std::endl;
      auto x_nrad_ = stan::math::eval(f(ret_val, y_nrad_, arena_dat, arena_dat_int, msgs));
      std::cout << "x_nrad: " << x_nrad_.val().transpose() << std::endl;
      //auto x_nrad_ = stan::math::eval(fy_(y_nrad_));
      x_nrad_.adj() = eta;
      std::cout << "y_adj1: " << arena_y.adj().transpose() << std::endl;
      stan::math::grad();
      std::cout << "y_adj2: " << arena_y.adj().transpose() << std::endl;
      std::cout << "y_nrad_: " << y_nrad_.adj().transpose() << std::endl;
      arena_y.adj() += y_nrad_.adj();
      std::cout << "y_adj3: " << arena_y.adj().transpose() << std::endl;
    }
  });

  return ret_type(ret);
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
    double relative_tolerance = 1e-10, double function_tolerance = 1e-6,
    long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
  return algebra_solver_powell(f, x, y, dat, dat_int, msgs, relative_tolerance,
                               function_tolerance, max_num_steps);
}

template <typename S, typename F, typename T,
          require_eigen_vector_t<T>* = nullptr>
Eigen::VectorXd algebra_solver_powell_(
    S& solver, const F& fx, const T& x, const Eigen::VectorXd& y,
    const std::vector<double>& dat, const std::vector<int>& dat_int,
    std::ostream* msgs, double relative_tolerance, double function_tolerance,
    long int max_num_steps) {  // NOLINT(runtime/int)
  const auto& x_eval = x.eval();
  const auto& x_val = (value_of(x_eval)).eval();

  // Compute theta_dbl
  Eigen::VectorXd theta_dbl = x_val;
  solver.parameters.xtol = relative_tolerance;
  solver.parameters.maxfev = max_num_steps;
  solver.solve(theta_dbl);

  // Check if the max number of steps has been exceeded
  if (solver.nfev >= max_num_steps) {
    throw_domain_error("algebra_solver", "maximum number of iterations",
                       max_num_steps, "(", ") was exceeded in the solve.");
  }

  // Check solution is a root
  double system_norm = fx.get_value(theta_dbl).stableNorm();
  if (system_norm > function_tolerance) {
    std::ostringstream message;
    message << "the norm of the algebraic function is " << system_norm
            << " but should be lower than the function "
            << "tolerance:";
    throw_domain_error("algebra_solver", message.str().c_str(),
                       function_tolerance, "",
                       ". Consider decreasing the relative tolerance and "
                       "increasing max_num_steps.");
  }

  return theta_dbl;
}

}  // namespace math
}  // namespace stan

#endif
