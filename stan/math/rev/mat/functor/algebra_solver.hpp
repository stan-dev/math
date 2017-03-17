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

    template <typename F, typename T1, typename T2>
    struct hybrj_functor_solver : stan::math::NLOFunctor<double> {
    private:
      F f;
      Eigen::MatrixXd J;

    public:
      hybrj_functor_solver(const F& f_param,
                           const Eigen::Matrix<T1, Eigen::Dynamic, 1> theta,
                           const Eigen::Matrix<T2, Eigen::Dynamic, 1> parms,
                           const std::vector<double> dat,
                           const std::vector<int> dat_int,
                           const std::string variable) :
                           f(theta, parms, dat, dat_int, variable) { }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
        stan::math::jacobian(f, x, fvec, J);
        return 0;
      }

      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = J;
        return 0;
      }

      template <typename T>
      inline 
      Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2, T>::type,
                     Eigen::Dynamic, 1>
        operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
        return f(x);
      }

      /*
      Eigen::MatrixXd get_jacobian(const Eigen::VectorXd &parms) {
        Eigen::VectorXd fvec;
        stan::math::jacobian(f, parms, fvec, J);
        return J;
      } */
    };

    namespace {
      template <typename F, typename T>
      class algebra_solver_vari_alloc : public chainable_alloc {
      private:
        inline void compute(const F& f,
                            const Eigen::Matrix<double,
                                                Eigen::Dynamic, 1>& x,
                            const Eigen::Matrix<double, 
                                                Eigen::Dynamic, 1>& parms,
                            const std::vector<double>& dat,
                            const std::vector<int>& dat_int) {
          hybrj_functor_solver<F, double, double>
            functor(f, x, parms, dat, dat_int, "theta");
          Eigen::HybridNonLinearSolver<hybrj_functor_solver<F, double, double> >
            solver(functor);
          Eigen::Matrix<double, Eigen::Dynamic, 1> theta_d = x;
          solver.solve(theta_d);

          for (int i = 0; i < theta_d.rows(); i++)
            theta_(i) = var(new vari(theta_d(i), false));
        }

      public:
        algebra_solver_vari_alloc(const F& f,
                                  const Eigen::Matrix<double,
                                    Eigen::Dynamic, 1>& x,
                                  const Eigen::Matrix<T,
                                    Eigen::Dynamic, 1>& parms,
                                  const std::vector<double>& dat,
                                  const std::vector<int>& dat_int)
          : f_(f), x_(x), parms_(parms), dat_(dat),
            dat_int_(dat_int), theta_(x_.rows()) {
          compute(f, x, value_of(parms), dat, dat_int);
        }

        F f_;
        Eigen::Matrix<double, Eigen::Dynamic, 1> x_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> parms_;
        std::vector<double> dat_;
        std::vector<int> dat_int_;
        Eigen::Matrix<T, Eigen::Dynamic, 1> theta_;
      };

      template <typename F, typename T>
      class algebra_solver_vari : public vari {
      protected:
        inline void chain(const F& f,
                          const Eigen::Matrix<double,
                            Eigen::Dynamic, 1>& x,
                          Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                          const std::vector<double>& dat,
                          const std::vector<int>& dat_int,
                          const Eigen::Matrix<double,
                            Eigen::Dynamic, 1>& adjTheta) {
          using Eigen::Matrix;
          using Eigen::Dynamic;
          using Eigen::VectorXd;
          using Eigen::MatrixXd;
          using Eigen::HybridNonLinearSolver;
          using std::vector;
          
          std::cout << "Call to algebra chain" << std::endl;  // TEST

          // Root of equation is stored in impl_->theta_ (need to acces val_)
          /*
          Matrix<double, Dynamic, 1> parms_val = value_of(parms);
          hybrj_functor_solver<F, double, double>
            functor(f, x, parms_val, dat, dat_int, "theta");
          HybridNonLinearSolver<hybrj_functor_solver<F, double, double> >
            solver(functor);
          Matrix<double, Dynamic, 1> theta_d = x;
          solver.solve(theta_d); */

          /*
          Eigen::MatrixXd Jf_x(theta_d.rows(), x.rows());
          Jf_x = functor.get_jacobian(theta_d);
          hybrj_functor_solver<F, double, double>
            functor_parm(f, theta_d, parms_val, dat, dat_int, "parms");
          Eigen::MatrixXd Jf_p(theta_d.rows(), parms_val.rows());  
          Jf_p = functor_parm.get_jacobian(parms_val);
          Eigen::MatrixXd Jx_p(x.rows(), parms_val.rows());
          Jx_p = - stan::math::mdivide_left(Jf_x, Jf_p);
          std::cout << Jx_p << std::endl; */

          int M = parms.rows();  // number of parameters
          int N = x.rows();  // number of states

          MatrixXd J_fy(N, M);
          MatrixXd J_fx(N, N);

          try {
            start_nested();

            vector<var> y_vars(M);
            // y_vars.reserve(M);
            for (int i = 0; i < M; i++) y_vars[i] = value_of(parms(i));        
            VectorXd y_input_dbl(M);
            for (int i = 0; i < M; i++) y_input_dbl(i) = y_vars[i].val();
            Matrix<var, Dynamic, 1> y_input_vars(M);
            for (int i = 0; i < M; i++) y_input_vars(i) = y_vars[i];
            vector<double> y_grad(M);

            vector<var> x_vars(N);
            // x_vars.reserve(N);
            for (int i = 0; i < N; i++)
              x_vars[i] = value_of(impl_->theta_(i));
            VectorXd x_input_dbl(N);
            for (int i = 0; i < N; i++) x_input_dbl(i) = x_vars[i].val();
            Matrix<var, Dynamic, 1> x_input_vars(N);
            for (int i = 0; i < N; i++) x_input_vars(i) = x_vars[i];
            vector<double> x_grad(N);

            hybrj_functor_solver<F, double, double>
              f_y(f, x_input_dbl, y_input_dbl, dat, dat_int, "parms");
            Matrix<var, Dynamic, 1> z_y = f_y(y_input_vars);

            hybrj_functor_solver<F, double, double>
              f_x(f, x_input_dbl, y_input_dbl, dat, dat_int, "theta");
            Matrix<var, Dynamic, 1> z_x = f_x(x_input_vars);

            for (int i = 0; i < N; i++) {
              set_zero_all_adjoints_nested();

              z_y(i).grad(y_vars, y_grad);
              for (int j = 0; j < M; j++) J_fy(i, j) = y_grad[j];

              z_x(i).grad(x_vars, x_grad);
              for (int j = 0; j < N; j++) J_fx(i, j) = x_grad[j];
            }
          } catch (const std::exception& e) {
            recover_memory_nested();
            throw;
          }

          recover_memory_nested();
          
          std::cout << "J_fy: " << std::endl << J_fy << std::endl;
          std::cout << "J_fx: " << std::endl << J_fx << std::endl;

          Eigen::MatrixXd J_xy(N, M);
          J_xy = - stan::math::mdivide_left(J_fx, J_fy);

          // Eigen::MatrixXd Jx_p(x.rows(), parms_val.rows());
          J_xy << 4, 5, 0,
                  0, 0, 1;

          std::cout << "algebra chain 1" << std::endl;

          for (int i = 0; i < adjTheta.rows(); i++)
            for (int j = 0; j < parms.rows(); j++)
              parms(j).vi_->adj_ += adjTheta(i) * J_xy(i, j);

          std::cout << "End of algebra chain" << std::endl;  // TEST
        }

      public:
        algebra_solver_vari(const F& f,
                            const Eigen::Matrix<double,
                              Eigen::Dynamic, 1>& x,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                            const std::vector<double>& dat,
                            const std::vector<int>& dat_int) : vari(0.0) {
          impl_
            = new algebra_solver_vari_alloc<F, T>(f, x, parms, dat, dat_int);
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

            chain(impl_->f_, impl_->x_, impl_->parms_, impl_->dat_,
                  impl_->dat_int_, adjTheta);
          }

          algebra_solver_vari_alloc<F, T> *impl_;
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
    template <typename T, typename F>
    inline Eigen::Matrix<var, Eigen::Dynamic, 1>
    algebra_solver(const F& f,
                   const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                   const Eigen::Matrix<T, Eigen::Dynamic, 1>& parms,
                   const std::vector<double>& dat,
                   const std::vector<int>& dat_int) {
      algebra_solver_vari<F, T> *baseVari
        = new algebra_solver_vari<F, T>(f, x, parms, dat, dat_int);

      return baseVari->impl_->theta_;
    }
  }
}
#endif
