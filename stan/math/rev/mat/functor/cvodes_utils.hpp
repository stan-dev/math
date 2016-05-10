#ifndef STAN_MATH_REV_MAT_FUNCTOR_CVODES_UTILS_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CVODES_UTILS_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/rev/mat/functor/ode_system.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

#include <vector>
#include <algorithm>
#include <string>

namespace stan {
  namespace math {

    // Noop error handler to silence CVodes error output
    extern "C"
    void cvodes_silent_err_handler(int error_code, const char *module,
                                   const char *function, char *msg,
                                   void *eh_data) {
    }

    void cvodes_check_flag(int flag, const std::string& func_name) {
      if (flag < 0) {
        std::ostringstream ss;
        ss << func_name << " failed with error flag " << flag;
        throw std::runtime_error(ss.str());
      }
    }

    void cvodes_set_options(void* cvodes_mem,
                            double rel_tol,
                            double abs_tol,
                            long int max_num_steps  // NOLINT(runtime/int)
                            ) {
      // Forward CVode errors to noop error handler
      CVodeSetErrHandlerFn(cvodes_mem, cvodes_silent_err_handler, 0);

      // Initialize solver parameters
      cvodes_check_flag(CVodeSStolerances(cvodes_mem, rel_tol, abs_tol),
                        "CVodeSStolerances");

      cvodes_check_flag(CVodeSetMaxNumSteps(cvodes_mem, max_num_steps),
                        "CVodeSetMaxNumSteps");

      double init_step = 0;
      cvodes_check_flag(CVodeSetInitStep(cvodes_mem, init_step),
                        "CVodeSetInitStep");

      long int max_err_test_fails = 20;  // NOLINT(runtime/int)
      cvodes_check_flag(CVodeSetMaxErrTestFails(cvodes_mem, max_err_test_fails),
                        "CVodeSetMaxErrTestFails");

      long int max_conv_fails = 50;  // NOLINT(runtime/int)
      cvodes_check_flag(CVodeSetMaxConvFails(cvodes_mem, max_conv_fails),
                        "CVodeSetMaxConvFails");
    }

    template<typename F, typename T_initial, typename T_param>
    class cvodes_ode_data {
      const std::vector<T_initial>& y0_;
      const std::vector<T_param>& theta_;
      const size_t N_;
      const size_t M_;
      const size_t param_var_ind_;

      const ode_system<F> ode_system_;

      typedef cvodes_ode_data<F, T_initial, T_param> ode_data;
      typedef stan::is_var<T_initial> initial_var;

    public:
      cvodes_ode_data(const std::vector<T_initial>& y0,
                      const F& f,
                      const std::vector<T_param>& theta,
                      const std::vector<double>& x,
                      const std::vector<int>& x_int,
                      std::ostream* msgs)
        : y0_(y0),
          theta_(theta),
          N_(y0.size()),
          M_(theta.size()),
          param_var_ind_(initial_var::value ? N_ : 0),
          ode_system_(f, stan::math::value_of(theta), x, x_int, msgs) {}

      // Static wrapper for CVode callback
      static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* user_data) {
        const ode_data* explicit_ode
          = reinterpret_cast<const ode_data*>(user_data);
        explicit_ode->rhs(NV_DATA_S(y), NV_DATA_S(ydot), t);
        return 0;
      }

      // Static wrapper for CVode callback
      static int ode_rhs_sens(int Ns, realtype t,
                              N_Vector y, N_Vector ydot,
                              N_Vector *yS, N_Vector *ySdot, void *user_data,
                              N_Vector tmp1, N_Vector tmp2) {
        const ode_data* explicit_ode
          = reinterpret_cast<const ode_data*>(user_data);

        const std::vector<double> y_vec(NV_DATA_S(y),
                                        NV_DATA_S(y) + explicit_ode->N_);

        explicit_ode->rhs_sens(explicit_ode->y0_,
                               explicit_ode->theta_,
                               t, y_vec, yS, ySdot);
        return 0;
      }

      // Static wrapper for CVode callback
      static int dense_jacobian(long int N,  // NOLINT(runtime/int)
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *user_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        const ode_data* explicit_ode
          = reinterpret_cast<const ode_data*>(user_data);
        return explicit_ode->dense_jacobian(NV_DATA_S(y), J, t);
      }

    private:
      void rhs(const double y[], double dy_dt[], double t) const {
        const std::vector<double> y_vec(y, y + N_);

        std::vector<double> dy_dt_vec(N_);
        ode_system_(t, y_vec, dy_dt_vec);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      int dense_jacobian(const double* y, DlsMat J, double t) const {
        const std::vector<double> y_vec(y, y + N_);

        Eigen::VectorXd fy(N_);
        // Eigen and CVODES use column major addressing
        Eigen::Map<Eigen::MatrixXd> Jy_map(J->data, N_, N_);

        ode_system_.jacobian(t, y_vec, fy, Jy_map);

        return 0;
      }
          inline void
      rhs_sens_initial(const Eigen::MatrixXd& Jy,
                       N_Vector *yS, N_Vector *ySdot) const {
          for (size_t m = 0; m < N_; m++) {
            // map NV_Vector to Eigen facilities
            Eigen::Map<Eigen::VectorXd> yS_eig(NV_DATA_S(yS[m]), N_);
            Eigen::Map<Eigen::VectorXd> ySdot_eig(NV_DATA_S(ySdot[m]), N_);

            ySdot_eig = Jy * yS_eig;
          }
      }

      inline void
      rhs_sens_param(const Eigen::MatrixXd& Jy,
                     const Eigen::MatrixXd& Jtheta,
                     N_Vector *yS, N_Vector *ySdot) const {
        for (size_t m = 0; m < M_; m++) {
          // map NV_Vector to Eigen facilities
          Eigen::Map<Eigen::VectorXd>
            yS_eig(NV_DATA_S(yS[param_var_ind_ + m]), N_);
          Eigen::Map<Eigen::VectorXd>
            ySdot_eig(NV_DATA_S(ySdot[param_var_ind_ + m]), N_);

          ySdot_eig = Jy * yS_eig + Jtheta.col(m);
        }
      }

      /**
       * Calculate the sensitivity RHS for varying initials and parameters.
       *
       * @param[in] initial var vector
       * @param[in] param var vector
       * @param[in] t time
       * @param[in] y state of the base ODE system
       * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
       * RHS
       * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
       */
      void rhs_sens(const std::vector<stan::math::var>& initial,
                    const std::vector<stan::math::var>& param,
                    const double t, const std::vector<double>& y,
                    N_Vector *yS, N_Vector *ySdot) const {
        Eigen::VectorXd dy_dt(N_);
        Eigen::MatrixXd Jy(N_, N_);
        Eigen::MatrixXd Jtheta(N_, M_);
        ode_system_.jacobian(t, y, dy_dt, Jy, Jtheta);
        rhs_sens_initial(Jy, yS, ySdot);
        rhs_sens_param(Jy, Jtheta, yS, ySdot);
      }

      /**
       * Calculate the sensitivity RHS for fixed initials and varying parameters.
       *
       * @param[in] initial double vector
       * @param[in] param var vector
       * @param[in] t time
       * @param[in] y state of the base ODE system
       * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
       * RHS
       * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
       */
      void rhs_sens(const std::vector<double>& initial,
                    const std::vector<stan::math::var>& param,
                    const double t, const std::vector<double>& y,
                    N_Vector *yS, N_Vector *ySdot) const {
        Eigen::VectorXd dy_dt(N_);
        Eigen::MatrixXd Jy(N_, N_);
        Eigen::MatrixXd Jtheta(N_, M_);
        ode_system_.jacobian(t, y, dy_dt, Jy, Jtheta);
        rhs_sens_param(Jy, Jtheta, yS, ySdot);
      }

      /**
       * Calculate the sensitivity RHS for varying initials and fixed parameters.
       *
       * @param[in] initial var vector
       * @param[in] param double vector
       * @param[in] t time
       * @param[in] y state of the base ODE system
       * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
       * RHS
       * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
       */
      void rhs_sens(const std::vector<stan::math::var>& initial,
                    const std::vector<double>& param,
                    const double t, const std::vector<double>& y,
                    N_Vector *yS, N_Vector *ySdot) const {
        Eigen::VectorXd dy_dt(N_);
        Eigen::MatrixXd Jy(N_, N_);
        ode_system_.jacobian(t, y, dy_dt, Jy);
        rhs_sens_initial(Jy, yS, ySdot);
      }

      /**
       * Calculate the empty sensitivity RHS.
       *
       * @param[in] initial double vector
       * @param[in] param double vector
       * @param[in] t time
       * @param[in] y state of the base ODE system
       * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
       * RHS
       * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
       */
      void rhs_sens(const std::vector<double>& initial,
                    const std::vector<double>& param,
                    const double t, const std::vector<double>& y,
                    N_Vector *yS, N_Vector *ySdot) const {
      }
    };
  }  // math
}  // stan

#endif
