#ifndef STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP
#define STAN_MATH_PRIM_MAT_FUN_ALGEBRA_SOLVER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include<stan/math/prim/mat/fun/dogleg.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/prim/mat/fun/inverse.hpp>

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

      int operator()(const Eigen::VectorXd &x, 
                     Eigen::VectorXd &fvec) {
        stan::math::jacobian(f1, x, fvec, J);
        return 0;
      }

      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = J;  //f2(x);
        return 0;
      }

      Eigen::MatrixXd get_jacobian(const Eigen::VectorXd &parms) {
        Eigen::VectorXd fvec;
        stan::math::jacobian(f1, parms, fvec, J);
        return J;
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
     template <typename F1, typename T1, typename T2>
     inline
     Eigen::Matrix<typename stan::return_type<T1, T2>::type,
       Eigen::Dynamic, Eigen::Dynamic>
     algebra_solver(const F1& f1,
                    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
                    const std::vector<double>& dat,
                    const std::vector<int>& dat_int) {
       typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;

       stan::math::hybrj_functor_solver<F1, T1, T2>
         functor(f1, x, parms, dat, dat_int, "theta");
       Eigen::HybridNonLinearSolver<hybrj_functor_solver<F1, T1, T2> >
         solver(functor);
       Eigen::Matrix<scalar, Eigen::Dynamic, 1> theta = x;
       solver.solve(theta);

       // CHECK - should the argument be x or theta ?
       // Do we care about the derivative w.r.t x.
       Eigen::MatrixXd Jf_x = functor.get_jacobian(theta);

       stan::math::hybrj_functor_solver<F1, T1, T2>
         functor_parm(f1, theta, parms, dat, dat_int, "parms");
       Eigen::MatrixXd Jf_p = functor_parm.get_jacobian(parms);
       
       Eigen::MatrixXd Jx_p = - stan::math::inverse(Jf_x) * Jf_p;

       // print statements
       /*
       using std::cout;
       using std::endl;
       cout << "Jf_x" << endl << Jf_x << endl;
       cout << "Jf_p" << endl << Jf_p << endl;
       cout << "inv(Jf_x)" << endl << stan::math::inverse(Jf_x) << endl;
       cout << "Jx_p" << endl << Jx_p << endl << endl;
       */
       
       return theta;
    }
  }
}

#endif