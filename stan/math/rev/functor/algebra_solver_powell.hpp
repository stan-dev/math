#ifndef STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP
#define STAN_MATH_REV_FUNCTOR_ALGEBRA_SOLVER_POWELL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
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
 * The vari class for the algebraic solver. We compute the  Jacobian of
 * the solutions with respect to the parameters using the implicit
 * function theorem. The call to Jacobian() occurs outside the call to
 * chain() -- this prevents malloc issues.
 */
template <typename Fs, typename F, typename T, typename Fx>
struct algebra_solver_vari : public vari {
  /** vector of parameters */
  vari** y_;
  /** number of parameters */
  int y_size_;
  /** number of unknowns */
  int x_size_;
  /** vector of solution */
  vari** theta_;
  /** Jacobian of the solution w.r.t parameters */
  double* Jx_y_;

  algebra_solver_vari(const Fs& fs, const F& f, const Eigen::VectorXd& x,
                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                      const std::vector<double>& dat,
                      const std::vector<int>& dat_int,
                      const Eigen::VectorXd& theta_dbl, Fx& fx,
                      std::ostream* msgs)
      : vari(theta_dbl(0)),
        y_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(y.size())),
        y_size_(y.size()),
        x_size_(x.size()),
        theta_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(x_size_)),
        Jx_y_(ChainableStack::instance_->memalloc_.alloc_array<double>(
            x_size_ * y_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;
    for (int i = 0; i < y.size(); ++i) {
      y_[i] = y(i).vi_;
    }

    theta_[0] = this;
    for (int i = 1; i < x.size(); ++i) {
      theta_[i] = new vari(theta_dbl(i), false);
    }

    // Compute the Jacobian and store in array, using the
    // implicit function theorem, i.e. Jx_y = Jf_y / Jf_x
    using f_y = hybrj_functor_solver<Fs, F, double, double>;
    Map<MatrixXd>(&Jx_y_[0], x_size_, y_size_)
        = -mdivide_left(fx.get_jacobian(theta_dbl),
                        f_y(fs, f, theta_dbl, value_of(y), dat, dat_int, msgs)
                            .get_jacobian(value_of(y)));
  }

  void chain() {
    for (int j = 0; j < y_size_; j++) {
      for (int i = 0; i < x_size_; i++) {
        y_[j]->adj_ += theta_[i]->adj_ * Jx_y_[j * x_size_ + i];
      }
    }
  }
};

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
  const auto& x_eval = x.eval();
  const auto& x_val = (value_of(x_eval)).eval();
  algebra_solver_check(x_val, y, dat, dat_int, function_tolerance,
                       max_num_steps);
  check_nonnegative("alegbra_solver", "relative_tolerance", relative_tolerance);
  // if (relative_tolerance < 0)
  //   invalid_argument("algebra_solver", "relative_tolerance,",
  //                    function_tolerance, "",
  //                    ", must be greater than or equal to 0");

  // Create functor for algebraic system
  using Fs = system_functor<F, double, double, true>;
  using Fx = hybrj_functor_solver<Fs, F, double, double>;
  Fx fx(Fs(), f, x_val, y, dat, dat_int, msgs);
  Eigen::HybridNonLinearSolver<Fx> solver(fx);

  // Check dimension unknowns equals dimension of system output
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       fx.get_value(x_val), "the vector of unknowns, x,", x);

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
  const auto& y_eval = y.eval();
  const auto& x_val = (value_of(x_eval)).eval();
  const auto& y_val = (value_of(y_eval)).eval();
  Eigen::VectorXd theta_dbl = algebra_solver_powell(
      f, x_eval, y_val, dat, dat_int, 0, relative_tolerance, function_tolerance,
      max_num_steps);

  using Fy = system_functor<F, double, double, false>;

  // TODO(charlesm93): a similar object gets constructed inside
  // the call to algebra_solver. Cache the previous result
  // and use it here (if possible).
  using Fs = system_functor<F, double, double, true>;
  using Fx = hybrj_functor_solver<Fs, F, double, double>;
  Fx fx(Fs(), f, x_val, y_val, dat, dat_int, msgs);

  // Construct vari
  algebra_solver_vari<Fy, F, value_type_t<T2>, Fx>* vi0
      = new algebra_solver_vari<Fy, F, value_type_t<T2>, Fx>(
          Fy(), f, x_val, y_eval, dat, dat_int, theta_dbl, fx, msgs);
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(x.size());
  theta(0) = var(vi0->theta_[0]);
  for (int i = 1; i < x.size(); ++i) {
    theta(i) = var(vi0->theta_[i]);
  }

  return theta;
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

}  // namespace math
}  // namespace stan

#endif
