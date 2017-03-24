#ifndef STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <stdexcept>  // TEST
#include <iostream>

namespace stan {
  namespace math {

    template <typename F, typename T1, typename T2>
    struct hybrj_functor_solver : stan::math::NLOFunctor<double> {
    private:
      F f_;
      Eigen::MatrixXd J_;

    public:
      hybrj_functor_solver(const F& f,
                           const Eigen::Matrix<T1, Eigen::Dynamic, 1> x,
                           const Eigen::Matrix<T2, Eigen::Dynamic, 1> y,
                           const std::vector<double> dat,
                           const std::vector<int> dat_int,
                           const bool x_is_dv) :
                           f_(x, y, dat, dat_int, x_is_dv) { }

      int operator()(const Eigen::VectorXd &dv, Eigen::VectorXd &fvec) {
        stan::math::jacobian(f_, dv, fvec, J_);
        return 0;
      }

      int df(const Eigen::VectorXd &dv, Eigen::MatrixXd &fjac) {
        fjac = J_;
        return 0;
      }

      Eigen::MatrixXd get_jacobian(const Eigen::VectorXd &y) {
        Eigen::VectorXd fvec;
        stan::math::jacobian(f_, y, fvec, J_);
        return J_;
      }
    };

    template <typename F, typename T, typename FX>
    struct algebra_solver_vari : public vari {
      vari** y_;
      int y_size_;
      int x_size_;
      vari** theta_;
      Eigen::MatrixXd Jx_y_;

      algebra_solver_vari(const F& f,
                          const Eigen::VectorXd x,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1> y,
                          const std::vector<double> dat,
                          const std::vector<int> dat_int,
                          const Eigen::VectorXd theta_dbl,
                          FX& fx)
        : vari(theta_dbl(0)),  // CHECK - calls chain()
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

        hybrj_functor_solver<F, double, double>
          fy(f, x, value_of(y), dat, dat_int, false);
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
     * which get passed into the algebraic system.
     *
     * @tparam F1 type of equation system function.
     * @tparam T type of scalars for parms.
     * @param[in] F1 Functor that evaluates the system of equations.
     * @param[in] x Vector of starting values.
     * @param[in] parms parameter vector for the equation system.
     * @param[in] dat continuous data vector for the equation system.
     * @param[in] dat_int integer data vector for the equation system.
     * @return Vector that solves the system of equations.
     */
    template <typename F, typename T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    algebra_solver(const F& f,
                   const Eigen::VectorXd& x,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int) {
      // Compute theta_dbl
      typedef hybrj_functor_solver<F, double, double> FX;
      FX fx(f, x, value_of(y), dat, dat_int, true);
      Eigen::HybridNonLinearSolver<FX> solver(fx);
      Eigen::VectorXd theta_dbl = x;
      solver.solve(theta_dbl);

      // Construct vari
      algebra_solver_vari<F, T, FX>* vi0
        = new algebra_solver_vari<F, T, FX>(f, x, y, dat, dat_int, theta_dbl, fx);
      Eigen::Matrix<T, Eigen::Dynamic, 1> theta(x.size());
      theta(0) = var(vi0);
      for (int i = 1; i < x.size(); ++i)
        theta(i) = var(vi0->theta_[i]);
      return theta;
     }
  }
}

#endif
