#ifndef STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP
#define STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>

namespace stan {
  namespace math {

    template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct NLOFunctor {
      typedef _Scalar Scalar;
      enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
      };

      typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
      typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
      typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

      const int m_inputs, m_values;

      NLOFunctor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
      NLOFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

      int inputs() const { return m_inputs; }
      int values() const { return m_values; }
    };

    template <typename F1>
    struct hybrj_functor : NLOFunctor<double> {
    private:
      F1 f1;
      Eigen::MatrixXd J;

    public:
      hybrj_functor(const F1& f1_param) {
        f1 = f1_param;
      }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
        stan::math::jacobian(f1, x, fvec, J);
        return 0;
      }
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = J; // FIXME: check that df() is never called before ()
        return 0;
      }
    };

    /**
     * @param x Vector of starting values
     * @param F1 Functor that evaluates the system of equations
     * @return Vector that solves the system of equations
     */
    template <typename F1>
    inline
    Eigen::VectorXd
    dogleg(const Eigen::Matrix<var, Eigen::Dynamic, 1>& x, const F1 f1) {

      stan::math::hybrj_functor<F1> functor(f1);
      Eigen::HybridNonLinearSolver<hybrj_functor<F1> > solver(functor);
      Eigen::VectorXd theta = stan::math::value_of(x);
      solver.solve(theta);
      return theta;
    }

  }
}
#endif
