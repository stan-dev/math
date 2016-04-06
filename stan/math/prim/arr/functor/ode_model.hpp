#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_ODE_MODEL_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_ODE_MODEL_HPP

#include <stan/math/rev/core.hpp>
#include <iostream>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Internal representation of ODE model object which provides
     * convenient jacobian functions to obtain gradients wrt to states
     * (S) and parameters (P). Can be used to provide analytic
     * Jacobians via partial template specialisation.
     */

    template<typename F>
    struct ode_model {
      const F& f_;
      const std::vector<double>& theta_;
      const std::vector<double>& x_;
      const std::vector<int>& x_int_;
      std::ostream* msgs_;

      ode_model(const F& f,
                const std::vector<double>& theta,
                const std::vector<double>& x,
                const std::vector<int>& x_int,
                std::ostream* msgs)
        : f_(f),
          theta_(theta),
          x_(x),
          x_int_(x_int),
          msgs_(msgs)
      {}

      void operator()(const std::vector<double>& y,
                      std::vector<double>& dy_dt,
                      const double t) const {
        dy_dt = f_(t, y, theta_, x_, x_int_, msgs_);
      }

      template <typename Derived1, typename Derived2>
      void
      jacobian_S(const double t,
                 const std::vector<double>& y,
                 Eigen::MatrixBase<Derived1>& fy,
                 Eigen::MatrixBase<Derived2>& Jy) const {
        // note: the function expects fy and Jy to have correct sizes!
        using Eigen::Matrix;
        using stan::math::var;
        using std::vector;
        stan::math::start_nested();
        try {
          vector<var> y_var(y.size());
          for (size_t k = 0; k < y.size(); ++k)
            y_var[k] = y[k];
          vector<var> fy_var = f_(t, y_var, theta_, x_, x_int_, msgs_);
          for (size_t i = 0; i < fy_var.size(); ++i)
            fy(i) = fy_var[i].val();
          for (size_t i = 0; i < fy_var.size(); ++i) {
            if (i > 0)
              stan::math::set_zero_all_adjoints_nested();
            grad(fy_var[i].vi_);
            for (size_t k = 0; k < y.size(); ++k)
              Jy(i, k) = y_var[k].adj();
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      template <typename Derived1, typename Derived2, typename Derived3>
      void
      jacobian_SP(const double t,
                  const std::vector<double>& y,
                  Eigen::MatrixBase<Derived1>& fy,
                  Eigen::MatrixBase<Derived2>& Jy,
                  Eigen::MatrixBase<Derived3>& Jtheta) const {
        using Eigen::Matrix;
        using Eigen::Dynamic;
        using stan::math::var;
        using std::vector;
        stan::math::start_nested();
        try {
          vector<var> y_var(y.size());
          for (size_t k = 0; k < y.size(); ++k)
            y_var[k] = y[k];
          vector<var> theta_var(theta_.size());
          for (size_t k = 0; k < theta_.size(); ++k)
            theta_var[k] = theta_[k];
          vector<var> z_var;
          z_var.reserve(y.size() + theta_.size());
          z_var.insert(z_var.end(),     y_var.begin(),     y_var.end());
          z_var.insert(z_var.end(), theta_var.begin(), theta_var.end());
          vector<var> fy_var = f_(t, y_var, theta_var, x_, x_int_, msgs_);
          for (size_t i = 0; i < fy_var.size(); ++i)
            fy(i) = fy_var[i].val();
          for (size_t i = 0; i < fy_var.size(); ++i) {
            if (i > 0)
              stan::math::set_zero_all_adjoints_nested();
            grad(fy_var[i].vi_);
            for (size_t k = 0; k < y.size(); ++k)
              Jy(i, k) = y_var[k].adj();
            for (size_t k = 0; k < theta_.size(); ++k)
              Jtheta(i, k) = theta_var[k].adj();
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }
    };
  }
}

#endif
