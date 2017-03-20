#ifndef STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/mat/functor/jacobian_test.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
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
    
    Eigen::MatrixXd get_jacobian_test(const Eigen::VectorXd &y) {
      Eigen::VectorXd fvec;
      stan::math::jacobian_test(f_, y, fvec, J_);
      return J_;
    }
  };

  namespace {
    template <typename F, typename T>
    class algebra_solver_vari_alloc : public chainable_alloc {
    private:
      inline void compute(const F& f,
                          const Eigen::VectorXd& x,
                          const Eigen::VectorXd& y,
                          const std::vector<double>& dat,
                          const std::vector<int>& dat_int) {
        hybrj_functor_solver<F, double, double>
          functor(f, x, y, dat, dat_int, true);
        Eigen::HybridNonLinearSolver<hybrj_functor_solver<F, double, double> >
          solver(functor);
        Eigen::VectorXd theta_dbl = x;
        solver.solve(theta_dbl);

        for (int i = 0; i < theta_dbl.rows(); i++)
          theta_(i) = var(new vari(theta_dbl(i), false));  // CHECK: second argument
          // theta_(i) = var(new vari(theta_dbl(i), true));
      }

    public:
      algebra_solver_vari_alloc(const F& f,
                                const Eigen::VectorXd& x,
                                const Eigen::Matrix<T,
                                  Eigen::Dynamic, 1>& y,
                                const std::vector<double>& dat,
                                const std::vector<int>& dat_int)
        : f_(f), x_(x), y_(y), dat_(dat), dat_int_(dat_int),
          theta_(x_.rows()) {  // CHECK
        compute(f, x, value_of(y), dat, dat_int);
      }

      F f_;
      Eigen::VectorXd x_;
      Eigen::Matrix<T, Eigen::Dynamic, 1> y_;
      std::vector<double> dat_;
      std::vector<int> dat_int_;
      Eigen::Matrix<T, Eigen::Dynamic, 1> theta_;
    };

    template <typename F, typename T>
    class algebra_solver_vari : public vari {
    protected:
      inline void chain(const F& f,
                        const Eigen::VectorXd& x,
                        const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                        const std::vector<double>& dat,
                        const std::vector<int>& dat_int,
                        const Eigen::VectorXd adjTheta) {
        assert(x.rows() == adjTheta.rows());
        assert(x.rows() == impl_->theta_.rows());

        Eigen::VectorXd theta_dbl = value_of(impl_->theta_);
        Eigen::VectorXd y_dbl = value_of(y);

        hybrj_functor_solver<F, double, double>
          functor(f, x, y_dbl, dat, dat_int, true);
        // CHECK: should I specify the dimensions of Jf_x ?
        // CHECK: can I skip this step, since Jf_x gets evaluated in
        // compute()?
        Eigen::MatrixXd Jf_x = functor.get_jacobian(theta_dbl);

        hybrj_functor_solver<F, double, double>
          functor_parm(f, theta_dbl, y_dbl, dat, dat_int, false);
        Eigen::MatrixXd Jf_y(x.rows(), y.rows());
        Jf_y = functor_parm.get_jacobian_test(y_dbl);  // TEST
        // BUG: the line below causes an error !!
        // Eigen::MatrixXd Jf_y = functor_parm.get_jacobian(y_dbl);

        // Eigen::MatrixXd Jf_y(x.rows(), y.rows());
        Jf_y << -4, -5, 0, 0, 0, -1;  // TEST

        Eigen::MatrixXd Jx_y = - stan::math::mdivide_left(Jf_x, Jf_y);
        assert(Jx_y.rows() == adjTheta.rows());
        assert(Jx_y.cols() == y.rows());

        Jx_y << 4, 5, 0, 0, 0, 1;  // TEST

        // CHECK - do the indices match? Is adjTheta initialized? (yes, yes)
        for (int i = 0; i < adjTheta.rows(); i++)
          for (int j = 0; j < y.rows(); j++)
            y(j).vi_->adj_ += adjTheta(i) * Jx_y(i, j);
        
        // TEST
        std::cout << "y: " << y(2).vi_->val_ << ":" << y(2).vi_->adj_
                  << " " << y(1).vi_->val_ << ":" << y(1).vi_->adj_
                  << " " << y(0).vi_->val_ << ":" << y(0).vi_->adj_
                  << std::endl;
        
      }

    public:
      algebra_solver_vari(const F& f,
                          const Eigen::VectorXd& x,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                          const std::vector<double>& dat,
                          const std::vector<int>& dat_int) : vari(0.0) {
        impl_
          = new algebra_solver_vari_alloc<F, T>(f, x, y, dat, dat_int);
      }

      virtual void chain() {
        Eigen::VectorXd adjTheta(impl_->theta_.rows());
        for (int i = 0; i < adjTheta.rows(); i++)
          adjTheta(i) = impl_->theta_(i).vi_->adj_;

        chain(impl_->f_, impl_->x_, impl_->y_, impl_->dat_,
              impl_->dat_int_, adjTheta);
      }

      algebra_solver_vari_alloc<F, T> *impl_;
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
    template <typename T, typename F>
    inline Eigen::Matrix<T, Eigen::Dynamic, 1>
    algebra_solver(const F& f,
                   const Eigen::VectorXd& x,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int) {
      algebra_solver_vari<F, T> *baseVari
        = new algebra_solver_vari<F, T>(f, x, y, dat, dat_int);

      return baseVari->impl_->theta_;
    }
  }

  }
}

#endif
