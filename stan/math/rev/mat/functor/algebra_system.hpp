#ifndef STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SYSTEM_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ALGEBRA_SYSTEM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
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

      system_functor(const F f,
                     const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
                     const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                     const std::vector<double>& dat,
                     const std::vector<int>& dat_int,
                     std::ostream* msgs,
                     const bool& x_is_dv)
        : f_() {
        x_ = x;
        y_ = y;
        dat_ = dat;
        dat_int_ = dat_int;
        msgs_ = msgs;
        x_is_dv_ = x_is_dv;
      }

      template <typename T>
      inline
      Eigen::Matrix<T, Eigen::Dynamic, 1>
      operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
        if (x_is_dv_)
          return f_(x, y_, dat_, dat_int_, msgs_);
        else
          return f_(x_, x, dat_, dat_int_, msgs_);
      }
    };

    /**
     * This functor has the structure required to call Eigen's solver.
     */
    template <typename FS, typename F, typename T0, typename T1>
    struct hybrj_functor_solver : stan::math::NLOFunctor<double> {
    private:
      FS fs_;
      int x_size_;
      Eigen::MatrixXd J_;

    public:
      hybrj_functor_solver(const FS& fs,
                           const F& f,
                           const Eigen::Matrix<T0, Eigen::Dynamic, 1> x,
                           const Eigen::Matrix<T1, Eigen::Dynamic, 1> y,
                           const std::vector<double> dat,
                           const std::vector<int> dat_int,
                           std::ostream* msgs,
                           const bool x_is_dv)
        : fs_(f, x, y, dat, dat_int, msgs, x_is_dv),
          x_size_(x.size()) { }

      int operator()(const Eigen::VectorXd &dv, Eigen::VectorXd &fvec) {
        stan::math::jacobian(fs_, dv, fvec, J_);
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

      Eigen::VectorXd get_value(const Eigen::VectorXd dv) {
        return fs_(dv);
      }
    };

  }
}

#endif
