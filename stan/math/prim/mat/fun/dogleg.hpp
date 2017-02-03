#ifndef STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP
#define STAN_MATH_PRIM_MAT_FUN_DOGLEG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <unsupported/Eigen/NonLinearOptimization>

namespace stan {
  namespace math {
/*
    struct NLOFunctor {
      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;
    }; */

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
      
      // int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
      // int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;
    };

    template <typename F1, typename F2>
    struct hybrj_functor : NLOFunctor<double> {
    private:
      F1 f1;
      F2 f2;

    public:
      hybrj_functor(const F1& f1_param, const F2& f2_param) {
        f1 = f1_param;
        f2 = f2_param;
      }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
        fvec = f1(x);
        return 0;
      }
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = f2(x);
        return 0;
      }
    };

    /**
     * @param x Vector of starting values
     * @param F1 Functor that evaluates the system of equations
     * @param F2 Functor that evaluates the Jacobian matrix
     * @return Vector that solves the system of equations
     */
    // template <typename F1, typename F2>
    template <typename F1, typename F2>
    inline
    Eigen::VectorXd
    dogleg(const Eigen::VectorXd& x, const F1 f1, const F2 f2) {

      stan::math::hybrj_functor<F1, F2> functor(f1, f2);
      Eigen::HybridNonLinearSolver<NLOFunctor<double> > solver(functor);
      Eigen::VectorXd theta = x, fvec;
      functor(x, fvec);
      /*std::cout << x << std::endl;
      std::cout << fvec << std::endl;
      Eigen::MatrixXd fjac;
      functor.df(x, fjac);
      std::cout << fjac << std::endl; */
      solver.solve(theta);
      return theta;
    }

  }
}
#endif
