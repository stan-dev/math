#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SYSTEM_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SYSTEM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
  namespace math {

    /**
     * A functor that allows us to treat either x or y as
     * the independent variable. This allows us to
     * choose with respect to which variable the
     * jacobian gets computed.
     */
    template <typename F, typename T0, typename T1>
    struct system_functor {
    private:
      F f_;
      Eigen::Matrix<T0, Eigen::Dynamic, 1> x_;
      Eigen::Matrix<T1, Eigen::Dynamic, 1> y_;
      std::vector<double> dat_;
      std::vector<int> dat_int_;
      std::ostream* msgs_;
      bool x_is_dv_;

    public:
      system_functor() { }

      system_functor(const F& f,
                     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
                     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                     const std::vector<double>& dat,
                     const std::vector<int>& dat_int,
                     std::ostream* msgs,
                     const bool& x_is_dv)
        : f_(), x_(x), y_(y), dat_(dat), dat_int_(dat_int), msgs_(msgs),
        x_is_dv_(x_is_dv) { }

      template <typename T>
      inline
      Eigen::Matrix<T, Eigen::Dynamic, 1>
      operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
        if (x_is_dv_)
          return f_(x, y_, dat_, dat_int_, msgs_);
        else
          return f_(x_, x, dat_, dat_int_, msgs_);
      }
    };

    /**
     * A structure which gets passed to Eigen's dogleg
     * algebraic solver.
     */
    template <typename T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct NLOFunctor {
      typedef Eigen::Matrix<T, NX, 1> InputType;
      typedef Eigen::Matrix<T, NY, 1> ValueType;
      typedef Eigen::Matrix<T, NY, NX>
        JacobianType;

      const int m_inputs, m_values;

      NLOFunctor() : m_inputs(NX),
      m_values(NY) {}
      NLOFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

      int inputs() const { return m_inputs; }
      int values() const { return m_values; }
    };

    /**
    * A functor which stores an algebraic system and
    * its gradient, and contains the class functions
    * required for Eigen's dogleg algebraic solver.
    *
    * Members:
    * f1 functor which returns output of algebraic system.
    * f2 gradient of f1 with respect to the unknowns for
    * which we solve.
    */
    template <typename F1, typename F2>
    struct hybrj_functor : NLOFunctor<double> {
    private:
      F1 f1_;
      F2 f2_;

    public:
      hybrj_functor(const F1& f1,
                    const F2& f2)
        : f1_(f1), f2_(f2) { }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
        fvec = f1_(x);
        return 0;
      }
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) {
        fjac = f2_(x);
        return 0;
      }
    };

    /**
     * A functor with the rquired operators to call Eigen's 
     * algebraic solver.
     */
    template <typename FS, typename F, typename T0, typename T1>
    struct hybrj_functor_solver : NLOFunctor<double> {
    private:
      FS fs_;
      int x_size_;
      Eigen::MatrixXd J_;

    public:
      hybrj_functor_solver(const FS& fs,
                           const F& f,
                           const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
                           const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                           const std::vector<double>& dat,
                           const std::vector<int>& dat_int,
                           std::ostream* msgs,
                           const bool& x_is_dv)
        : fs_(f, x, y, dat, dat_int, msgs, x_is_dv),
          x_size_(x.size()) { }

      int operator()(const Eigen::VectorXd &dv, Eigen::VectorXd &fvec) {
        jacobian(fs_, dv, fvec, J_);
        return 0;
      }

      int df(const Eigen::VectorXd &dv, Eigen::MatrixXd &fjac) {
        fjac = J_;
        return 0;
      }

      Eigen::MatrixXd get_jacobian(const Eigen::VectorXd &dv) {
        Eigen::VectorXd fvec;
        stan::math::jacobian(fs_, dv, fvec, J_);
        return J_;
      }

      Eigen::VectorXd get_value(const Eigen::VectorXd &dv) {
        return fs_(dv);
      }
    };

  }
}

#endif
