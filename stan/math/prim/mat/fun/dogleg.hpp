#ifndef STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP
#define STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <unsupported/Eigen/NonLinearOptimization>

namespace stan {
  namespace math {

    struct NLOFunctor {
      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;
    };

    /**
     * @param x Vector of starting values
     * @param F1 Functor that evaluates the system of equations
     * @param F2 Functor that evaluates the Jacobian matrix
     * @return Vector that solves the system of equations
     */
    template <typename F1, typename F2>
    inline
    Eigen::VectorXd
    dogleg(const Eigen::VectorXd& x, const F1, const F2) {

      struct hybrj_functor : NLOFunctor {

        int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
          fvec = F1(x);
          return 0;
        }
        int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
          fjac = F2(x);
          return 0;
        }
      };

      hybrj_functor functor;
      Eigen::HybridNonLinearSolver<NLOFunctor> solver(functor);
      Eigen::VectorXd theta = x;
      solver.solve(theta);
      return theta;
    }
  }
}
#endif
