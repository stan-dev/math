#ifndef STAN_MATH_REV_MAT_FUNCTOR_ODE_MODEL_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ODE_MODEL_HPP

#include <stan/math/rev/core.hpp>
#include <iostream>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Internal representation of ODE model object which provides
     * convenient Jacobian functions to obtain gradients wrt to states
     * (S) and parameters (P). Can be used to provide analytic
     * Jacobians via partial template specialisation.
     *
     * @tparam F type of functor for the base ode system.
     */
    template<typename F>
    struct ode_model {
      const F& f_;
      const std::vector<double>& theta_;
      const std::vector<double>& x_;
      const std::vector<int>& x_int_;
      std::ostream* msgs_;

      /**
       * Construct an ODE model with the specified base ODE system,
       * parameters, data, and a message stream.
       *
       * @param[in] f the base ODE system functor.
       * @param[in] theta parameters of the ode.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in, out] msgs stream to which messages are printed.
       */
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

      /**
       * Calculate the RHS of the ODE
       *
       * @param[in] y state of the ode system at time t.
       * @param[out] dy_dt ODE RHS
       * @param[in]  t time.
       */
      inline
      void operator()(const std::vector<double>& y,
                      std::vector<double>& dy_dt,
                      const double t) const {
        dy_dt = f_(t, y, theta_, x_, x_int_, msgs_);
      }

      /**
       * Return the number of parameters
       *
       * @return size of the theta vector
       */
      size_t size_param() const {
        return theta_.size();
      }

      /**
       * Calculate the Jacobian of the ODE RHS wrt to states y. The
       * function expects the output objects to have correct sizes,
       * i.e. dy_dt must be length N and Jy a NxN matrix (N states, M
       * parameters).
       *
       * @param[in] t time.
       * @param[in] y state of the ode system at time t.
       * @param[out] dy_dt ODE RHS
       * @param[out] Jy Jacobian of ODE RHS wrt to y.
       */
      template <typename Derived1, typename Derived2>
      void
      jacobian_S(const double t,
                 const std::vector<double>& y,
                 Eigen::MatrixBase<Derived1>& dy_dt,
                 Eigen::MatrixBase<Derived2>& Jy) const {
        using Eigen::Matrix;
        using stan::math::var;
        using std::vector;
        vector<double> grad(y.size());
        Eigen::Map<Eigen::RowVectorXd> grad_eig(&grad[0], y.size());
        try {
          stan::math::start_nested();
          vector<var> y_var(y.begin(), y.end());
          vector<var> dy_dt_var = f_(t, y_var, theta_, x_, x_int_, msgs_);
          for (size_t i = 0; i < dy_dt_var.size(); ++i) {
            dy_dt(i) = dy_dt_var[i].val();
            stan::math::set_zero_all_adjoints_nested();

            dy_dt_var[i].grad(y_var, grad);
            Jy.row(i) = grad_eig;
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
      }

      /**
       * Calculate the Jacobian of the ODE RHS wrt to states y and
       * parameters theta. The function expects the output objects to
       * have correct sizes, i.e. dy_dt must be length N, Jy a NxN
       * matrix and Jtheta a NxM matrix (N states, M parameters).
       *
       * @param[in] t time.
       * @param[in] y state of the ode system at time t.
       * @param[out] dy_dt ODE RHS
       * @param[out] Jy Jacobian of ODE RHS wrt to y.
       * @param[out] Jtheta Jacobian of ODE RHS wrt to theta.
       */
      template <typename Derived1, typename Derived2, typename Derived3>
      void
      jacobian_SP(const double t,
                  const std::vector<double>& y,
                  Eigen::MatrixBase<Derived1>& dy_dt,
                  Eigen::MatrixBase<Derived2>& Jy,
                  Eigen::MatrixBase<Derived3>& Jtheta) const {
        using Eigen::Matrix;
        using Eigen::Dynamic;
        using stan::math::var;
        using std::vector;
        vector<double> grad(y.size() + theta_.size());
        Eigen::Map<Eigen::RowVectorXd> grad_eig(&grad[0],
                                                y.size() + theta_.size());
        try {
          stan::math::start_nested();
          vector<var> y_var(y.begin(), y.end());
          vector<var> theta_var(theta_.begin(), theta_.end());
          vector<var> z_var;
          z_var.reserve(y.size() + theta_.size());
          z_var.insert(z_var.end(),     y_var.begin(),     y_var.end());
          z_var.insert(z_var.end(), theta_var.begin(), theta_var.end());
          vector<var> dy_dt_var = f_(t, y_var, theta_var, x_, x_int_, msgs_);
          for (size_t i = 0; i < dy_dt_var.size(); ++i) {
            dy_dt(i) = dy_dt_var[i].val();
            stan::math::set_zero_all_adjoints_nested();
            dy_dt_var[i].grad(z_var, grad);

            Jy.row(i) = grad_eig.leftCols(y.size());
            Jtheta.row(i) = grad_eig.rightCols(theta_.size());
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
