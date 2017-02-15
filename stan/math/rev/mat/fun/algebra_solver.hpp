#ifndef STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include<stan/math/prim/mat/fun/dogleg.hpp>
#include <unsupported/Eigen/NonLinearOptimization>

namespace stan {
  namespace math {

    template <typename F1, typename F2, typename T>
    struct hybrj_functor_solver : stan::math::NLOFunctor<double> {
    private:
      F1 f1;
      F2 f2;
//      std::vector<T> parms;
//      std::vector<double> dat;
//      std::vector<int> dat_int;

    public:
      hybrj_functor_solver(const F1& f1_param,
                           const F2& f2_param,
                           const std::vector<T> parms_param,
                           const std::vector<double> dat_param,
                           const std::vector<int> dat_int_param) :
                           f1(parms_param, dat_param, dat_int_param),
                           f2(parms_param, dat_param, dat_int_param) {
//        f1 = f1_param(parms_param, dat_param, dat_int_param);
//        f2 = f2_param(parms_param, dat_param, dat_int_param);
//        parms = parms_param;
//        dat = dat_param;
//        dat_int = dat_int_param;
      }

      int operator()(const Eigen::VectorXd &x, 
                     Eigen::VectorXd &fvec) {
        fvec = f1(x);
        return 0;
      }
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = f2(x);
        return 0;
      }
    };
    
    /**
     * @tparam F1 type of equation system function.
     * @tparam F2 type of Jacobian system function.
     * @tparam T1 type of scalars for initial guess.
     * @tparam T2 type of scalars for parms.
     * @param[in] x Vector of starting values.
     * @param[in] F1 Functor that evaluates the system of equations.
     * @param[in] F2 Functor that evaluates the Jacobian matrix.
     * @param[in] parms parameter vector for the equation system.
     * @param[in] dat continuous data vector for the equation system.
     * @param[in] dat_int integer data vector for the equation system.
     * @return Vector that solves the system of equations.
     */
     template <typename F1, typename F2, typename T1, typename T2>
     inline
     Eigen::Matrix<typename stan::return_type<T1, T2>::type,
       Eigen::Dynamic, Eigen::Dynamic>
     algebra_solver(const F1& f1,
                    const F2& f2,
                    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                    const std::vector<T2>& parms,
                    const std::vector<double>& dat,
                    const std::vector<int>& dat_int) {
       stan::math::hybrj_functor_solver<F1, F2, T2>
         functor(f1, f2, parms, dat, dat_int);
       Eigen::HybridNonLinearSolver<hybrj_functor_solver<F1, F2, T2> >
         solver(functor);
       Eigen::Matrix<T1, Eigen::Dynamic, 1> theta = x;  // CHECK - element type for theta
       solver.solve(theta);
       return theta;
    }
  }
}

#endif