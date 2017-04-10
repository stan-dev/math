#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * The vari class for the algebraic solver. We compute the  Jacobian of
     * the solutions with respect to the parameters using the implicit
     * function theorem. The call to Jacobian() occurs outside the call to
     * chain() -- this prevents malloc issues.
     *
     * @member y_ the vector of parameters
     * @member y_size_ the number of parameters
     * @member x_size_ the number of unknowns
     * @member theta_ the vector of solution
     * @member Jx_y_ the Jacobian of the solution with respect to the parameters.
     */
    template <typename FS, typename F, typename T, typename FX>
    struct algebra_solver_vari : public vari {
      vari** y_;
      int y_size_;
      int x_size_;
      vari** theta_;
      Eigen::MatrixXd Jx_y_;

      algebra_solver_vari(const FS& fs,
                          const F& f,
                          const Eigen::VectorXd x,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1> y,
                          const std::vector<double> dat,
                          const std::vector<int> dat_int,
                          const Eigen::VectorXd theta_dbl,
                          FX& fx,
                          std::ostream* msgs)
        : vari(theta_dbl(0)),
          y_(ChainableStack::memalloc_.alloc_array<vari*>(y.size())),
          y_size_(y.size()),
          x_size_(x.size()),
          theta_(ChainableStack::memalloc_.alloc_array<vari*>(x.size())) {
        for (int i = 0; i < y.size(); ++i)
          y_[i] = y(i).vi_;

        theta_[0] = this;
        for (int i = 1; i < x.size(); ++i)
          theta_[i] = new vari(theta_dbl(i), false);

        // Compute the Jacobian
        Eigen::MatrixXd Jf_x = fx.get_jacobian(theta_dbl);
        hybrj_functor_solver<FS, F, double, double>
          fy(fs, f, theta_dbl, value_of(y), dat, dat_int, msgs, false);
        Eigen::MatrixXd Jf_y = fy.get_jacobian(value_of(y));

        Jx_y_ = - stan::math::mdivide_left(Jf_x, Jf_y);
      }

      void chain() {
        for (int i = 0; i < x_size_; i++)
          for (int j = 0; j < y_size_; j++)
            y_[j]->adj_ += theta_[i]->adj_ * Jx_y_(i, j);
      }
    };

    /**
     * Return the solutions for the specified system of algebraic
     * equations given an initial guess, and parameters and data,
     * which get passed into the algebraic system. The user can
     * also specify the relative tolerance (xtol in Eigen's code),
     * the absolute tolerance, and the maximum number of steps
     * (maxfev in Eigen's code).
     *
     * Throw an exception if the norm of f(x), where f is the
     * output of the algebraic system and x the proposed solution,
     * is greater than the absolute tolerance. We here use the
     * norm as a metric of how far we are from the 0.
     *
     * @tparam F1 type of equation system function.
     * @tparam T type of scalars for parms.
     * @param[in] F1 Functor that evaluates the system of equations.
     * @param[in] x Vector of starting values.
     * @param[in] y parameter vector for the equation system.
     * @param[in] dat continuous data vector for the equation system.
     * @param[in] dat_int integer data vector for the equation system.
     * @param[in] relative_tolerance determines the convergence criteria 
     *            for the solution.
     * @param[in] absolute_tolerance determines whether roots are acceptable.
     * @param[in] max_num_steps  maximum number of function evaluations.
     * @return theta Vector of solutions to the system of equations.
     */
    template <typename F, typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    algebra_solver(const F& f,
                   const Eigen::VectorXd& x,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int,
                   std::ostream* msgs = 0,
                   double relative_tolerance = 1e-10,
                   double absolute_tolerance = 1e-6,
                   long int max_num_steps = 1e+3) {  // NOLINT(runtime/int)
      Eigen::VectorXd theta_dbl = algebra_solver(f, x, value_of(y), dat,
                                                 dat_int, 0,
                                                 relative_tolerance,
                                                 absolute_tolerance,
                                                 max_num_steps);

      // FIX ME - the next three lines are redeundant (they occur in the prim
      // algebra solver), but shouldn't be too expensive.
      typedef system_functor<F, double, double> FS;
      typedef hybrj_functor_solver<FS, F, double, double> FX;
      FX fx(FS(), f, x, value_of(y), dat, dat_int, msgs, true);

      // Construct vari
      algebra_solver_vari<FS, F,  T, FX>* vi0
        = new algebra_solver_vari<FS, F, T, FX>(FS(), f, x, y, dat, dat_int,
                                                theta_dbl, fx, msgs);
      Eigen::Matrix<T, Eigen::Dynamic, 1> theta(x.size());
      theta(0) = var(vi0);
      for (int i = 1; i < x.size(); ++i)
      theta(i) = var(vi0->theta_[i]);

      return theta;
     }
  }
}

#endif
