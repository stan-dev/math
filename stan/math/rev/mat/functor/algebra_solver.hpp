#ifndef STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <stan/math/prim/mat/fun/mdivide_left.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <iostream>

namespace stan {
  namespace math {

    template <typename F1, typename T1, typename T2>
    struct hybrj_functor_solver : stan::math::NLOFunctor<double> {
    private:
      F1 f1;
      Eigen::MatrixXd J;

    public:
      hybrj_functor_solver(const F1& f1_param,
                           const Eigen::Matrix<T1, Eigen::Dynamic, 1> theta,
                           const Eigen::Matrix<T2, Eigen::Dynamic, 1> parms,
                           const std::vector<double> dat,
                           const std::vector<int> dat_int,
                           const std::string variable) :
                           f1(theta, parms, dat, dat_int, variable) { }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
        stan::math::jacobian(f1, x, fvec, J);
        return 0;
      }

      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = J;
        return 0;
      }

      Eigen::MatrixXd get_jacobian(const Eigen::VectorXd &parms) {
        Eigen::VectorXd fvec;
        stan::math::jacobian(f1, parms, fvec, J);
        return J;
      }
    };

    namespace {
      template <typename F1, typename T>
      class algebra_solver_vari_alloc : public chainable_alloc {
      private:
        inline void compute(const F1& f1,
                            const Eigen::Matrix<double,
                              Eigen::Dynamic, 1>& x,
                            const Eigen::Matrix<double,
                              Eigen::Dynamic, 1>& parms,
                            const std::vector<double>& dat,
                            const std::vector<int>& dat_int) {
          hybrj_functor_solver<F1, double, double>
            functor(f1, x, parms, dat, dat_int, "theta");
          Eigen::HybridNonLinearSolver<hybrj_functor_solver<F1, double, double> >
            solver(functor);
          Eigen::Matrix<double, Eigen::Dynamic, 1> theta_d = x;
          solver.solve(theta_d);

          for (int i = 0; i < theta_d.rows(); i++)
            theta_(i) = var(new vari(theta_d(i), false));
        }

      public:
        algebra_solver_vari_alloc(const F1& f1,
                                  const Eigen::Matrix<double,
                                    Eigen::Dynamic, 1>& x,
                                  const Eigen::Matrix<T,
                                    Eigen::Dynamic, 1>& parms,
                                  const std::vector<double>& dat,
                                  const std::vector<int>& dat_int)
          : f1_(f1), x_(x), parms_(parms), dat_(dat),
            dat_int_(dat_int), theta_(x_.rows()) {
          compute(f1, x, value_of(parms), dat, dat_int);
        }

        F1 f1_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> x_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> parms_;
        std::vector<double> dat_;
        std::vector<int> dat_int_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> theta_;
      };

      template <typename F1, typename T>
      class algebra_solver_vari : public vari {
      protected:
        inline void chain(const F1& f1,
                          const Eigen::Matrix<double,
                            Eigen::Dynamic, 1>& x,
                          Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                          const std::vector<double>& dat,
                          const std::vector<int>& dat_int,
                          const Eigen::Matrix<double,
                            Eigen::Dynamic, 1>& adjTheta) {
          std::cout << "Call to algebra chain" << std::endl;  // TEST

          Eigen::MatrixXd Jx_p;  // CHECK - do I need default values?

          // find roots of the equation
          Eigen::Matrix<double, Eigen::Dynamic, 1>
            parms_val = value_of(parms);
          hybrj_functor_solver<F1, double, double>
            functor(f1, x, parms_val, dat, dat_int, "theta");
          Eigen::HybridNonLinearSolver<hybrj_functor_solver<F1, double, double> >
            solver(functor);
          Eigen::Matrix<double, Eigen::Dynamic, 1> theta = x;
          solver.solve(theta);

          // compute Jacobian
          Eigen::MatrixXd Jf_x = functor.get_jacobian(theta);

          Eigen::VectorXd p_dummy(3);
          p_dummy << 1, 5, 3;
          Eigen::VectorXd x_dummy(2);
          x_dummy << 6, 8;

          hybrj_functor_solver<F1, double, double>
            functor_parm(f1, theta, parms_val, dat, dat_int, "parms");

          // hybrj_functor_solver<F1, double, double>
            // functor_test(f1, x_dummy, p_dummy, dat, dat_int, "theta");
            // functor_test(f1, x_dummy, p_dummy, dat, dat_int, "parms");

          // Eigen::MatrixXd Jf_test = functor_test.get_jacobian(x_dummy);
          // Eigen::MatrixXd Jf_test = functor_test.get_jacobian(p_dummy);

          Eigen::MatrixXd Jf_p = functor_parm.get_jacobian(parms_val);
          Jx_p = - stan::math::mdivide_left(Jf_x, Jf_p);

          for (int i = 0; i < adjTheta.rows(); i++)
            for (int j = 0; j < parms.rows(); j++)
              parms(j).vi_->adj_ += adjTheta(i) * Jx_p(i, j);

          std::cout << "End of algebra chain" << std::endl;  // TEST
        }

      public:
        algebra_solver_vari(const F1& f1,
                            const Eigen::Matrix<double,
                              Eigen::Dynamic, 1>& x,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                            const std::vector<double>& dat,
                            const std::vector<int>& dat_int) : vari(0.0) {
          impl_
            = new algebra_solver_vari_alloc<F1, T>(f1, x, parms, dat, dat_int);
          }

          virtual void chain() {
            Eigen::Matrix<double, Eigen::Dynamic, 1>
              adjTheta(impl_->theta_.rows());

            std::cout << "theta adjoint: ";

            for (int i = 0; i < impl_->theta_.rows(); i++) {
              adjTheta(i) = impl_->theta_(i).vi_->adj_;
              std::cout << adjTheta(i) << " ";
            }
            std::cout << std::endl;

            chain(impl_->f1_, impl_->x_, impl_->parms_, impl_->dat_,
                  impl_->dat_int_, adjTheta);
          }

          algebra_solver_vari_alloc<F1, T> *impl_;
      };
    }

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
    template <typename T, typename F1>
    inline Eigen::Matrix<var, Eigen::Dynamic, 1>
    algebra_solver(const F1& f1,
                   const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int) {
      algebra_solver_vari<F1, T> *baseVari
        = new algebra_solver_vari<F1, T>(f1, x, parms, dat, dat_int);

      return baseVari->impl_->theta_;
    }
  }
}
#endif
