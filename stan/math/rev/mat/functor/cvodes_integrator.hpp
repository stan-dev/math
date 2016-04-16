#ifndef STAN_MATH_REV_MAT_FUNCTOR_CVODES_INTEGRATOR_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_CVODES_INTEGRATOR_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/rev/mat/functor/ode_model.hpp>

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
    void silent_err_handler_cvodes(int error_code, const char *module,
                                   const char *function, char *msg,
                                   void *eh_data) {
    }

    /**
     * CVODES integrator which supports stiff and non-stiff
     * integration with optional support for sensitivites for the
     * initial values and/or the parameters.
     * 
     * @tparam F type of functor for the base ode system.
    */
    template <typename F>
    class cvodes_integrator {
    private:
      const ode_model<F>& ode_model_;
      const std::vector<double>& y0_dbl_;
      double t0_;
      const bool initial_var_;
      const bool param_var_;
      const size_t M_;
      const size_t N_;
      const size_t S_;
      const size_t size_;
      const size_t param_var_ind_;
      void* cvode_mem_;
      std::vector<double> state_;
      N_Vector  cvode_state_;
      N_Vector *cvode_state_sens_;

      typedef cvodes_integrator<F> ode;

      void check_flag(int flag, const std::string& func_name) const {
        if (flag < 0) {
          std::ostringstream ss;
          ss << func_name << " failed with error flag " << flag;
          throw std::runtime_error(ss.str());
        }
      }

      void set_cvode_options(double rel_tol,
                             double abs_tol,
                             long int max_num_steps  // NOLINT(runtime/int)
                             ) {
        // Forward CVode errors to noop error handler
        CVodeSetErrHandlerFn(cvode_mem_, silent_err_handler_cvodes, 0);

        check_flag(CVodeInit(cvode_mem_, &ode::ode_rhs, t0_, cvode_state_),
                    "CVodeInit");

        // Assign pointer to this as user data
        check_flag(CVodeSetUserData(cvode_mem_,
                                    reinterpret_cast<void*>(this)),
                   "CVodeSetUserData");

        // Initialize solver parameters
        check_flag(CVodeSStolerances(cvode_mem_, rel_tol, abs_tol),
                   "CVodeSStolerances");

        check_flag(CVodeSetMaxNumSteps(cvode_mem_, max_num_steps),
                   "CVodeSetMaxNumSteps");

        double init_step = 0;
        check_flag(CVodeSetInitStep(cvode_mem_, init_step),
                   "CVodeSetInitStep");

        long int max_err_test_fails = 20;  // NOLINT(runtime/int)
        check_flag(CVodeSetMaxErrTestFails(cvode_mem_, max_err_test_fails),
                    "CVodeSetMaxErrTestFails");

        long int max_conv_fails = 50;  // NOLINT(runtime/int)
        check_flag(CVodeSetMaxConvFails(cvode_mem_, max_conv_fails),
                   "CVodeSetMaxConvFails");
      }

    public:
      /**
       * Construct CVODES integrator for an ODE model with initial
       * state, initial time, flag for sensitivity calculation of the
       * initial values and for the parameters, integrator options
       * such as rel+abs tolerance and maximum number of steps. The
       * integrator supports non-stiff (Adams-Moulton) and stiff (BDF)
       * integration with optional stability detection.
       *
       * The integrator creates as output a vector of vector
       * format. The outer vector is over the time-points where the
       * solution is requested and the inner the states. The order of
       * the states with is always
       *
       * \f[
       * (y, \frac{\partial y}{\partial y_0}, \frac{\partial y}{\partial \theta}).
       * \f]
       *
       * While the first N states correspond to y, all the remaining
       * are optional and depend on flags given to the constructor.
       *
       * @param[in] ode_model functor.
       * @param[in] y0 initial state of the base ode.
       * @param[in] t0 initial time of the base ode.
       * @param[in] initial_var flag if sensitivities of initals are needed
       * @param[in] param_var flag if sensitivities of initals are needed
       * @param[in] rel_tol Relative tolerance of solver.
       * @param[in] abs_tol Absolute tolerance of solver.
       * @param[in] max_num_steps Maximum number of solver steps.
       * @param[in] solver used solver (0=non-stiff, 1=stiff, 2=stiff
       * with STALD)
       */
      cvodes_integrator(const ode_model<F>& ode_model,
                        const std::vector<double>& y0,
                        double t0,
                        bool initial_var,
                        bool param_var,
                        double rel_tol,
                        double abs_tol,
                        long int max_num_steps,  // NOLINT(runtime/int)
                        size_t solver)
        : ode_model_(ode_model),
          y0_dbl_(y0),
          t0_(t0),
          initial_var_(initial_var),
          param_var_(param_var),
          M_(ode_model.size_param()),
          N_(y0.size()),
          S_((initial_var ? N_ : 0) + (param_var ? M_ : 0)),
          size_(N_ * (S_+1)),
          param_var_ind_(initial_var ? N_ : 0),
          cvode_mem_(NULL),
          state_(y0),
          cvode_state_(N_VMake_Serial(N_,
                                      &state_[0])),
          cvode_state_sens_(NULL) {
        stan::math::check_bounded("cvodes_integrator",
                                  "solver", solver,
                                  static_cast<size_t>(0),
                                  static_cast<size_t>(2));

        // Instantiate CVode memory
        if (solver == 0) {
          cvode_mem_ = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

          if (cvode_mem_ == NULL)
            throw std::runtime_error("CVodeCreate failed to allocate memory");

          set_cvode_options(rel_tol, abs_tol, max_num_steps);

        } else if (solver == 1 || solver == 2) {
          cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);

          if (cvode_mem_ == NULL)
            throw std::runtime_error("CVodeCreate failed to allocate memory");

          set_cvode_options(rel_tol, abs_tol, max_num_steps);

          // enable stability limit detection for solver==2
          if (solver == 2)
            check_flag(CVodeSetStabLimDet(cvode_mem_, 1), "CVodeSetStabLimDet");

          // for the stiff solvers we need to reserve additional
          // memory and provide a Jacobian function call
          check_flag(CVDense(cvode_mem_, N_), "CVDense");
          check_flag(CVDlsSetDenseJacFn(cvode_mem_, &ode::dense_jacobian),
                      "CVDlsSetDenseJacFn");
        }

        // initialize forward sensitivity system of CVODES as needed
        if (S_ > 0) {
          cvode_state_sens_ = N_VCloneVectorArray_Serial(S_, cvode_state_);

          for (size_t s = 0; s < S_ ; s++)
            N_VConst(RCONST(0.0), cvode_state_sens_[s]);

          // for the case with varying initials, the first N_
          // sensitivity systems correspond to the initials which have
          // as initial the identity matrix
          if (initial_var_) {
            for (size_t n = 0; n < N_; n++) {
              NV_Ith_S(cvode_state_sens_[n], n) = 1.0;
            }
          }

          check_flag(CVodeSensInit(cvode_mem_, static_cast<int>(S_),
                                   CV_STAGGERED,
                                   &ode::ode_rhs_sens, cvode_state_sens_),
                     "CVodeSensInit");

          check_flag(CVodeSensEEtolerances(cvode_mem_),
                     "CVodeSensEEtolerances");
        }
      }

      ~cvodes_integrator() {
        // N_VDestroy_Serial should noop because
        // N_VMake_Serial sets own_data to false
        N_VDestroy_Serial(cvode_state_);
        if (cvode_state_sens_ != NULL)
          N_VDestroyVectorArray_Serial(cvode_state_sens_, S_);
        CVodeFree(&cvode_mem_);
      }

      /**
       * Forward the CVode array-based call from
       * the ODE RHS to the Stan vector-based call
       *
       * @param[in] y state of the base ODE system
       * @param[out] dy_dt ODE RHS at time t
       * @param[in] t time
       */
      void rhs(const double y[], double dy_dt[], double t) const {
        const std::vector<double> y_vec(y, y + N_);

        std::vector<double> dy_dt_vec;
        ode_model_(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* f_data) {
        ode const* explicit_ode = reinterpret_cast<ode const*>(f_data);
        explicit_ode->rhs(NV_DATA_S(y), NV_DATA_S(ydot), t);
        return 0;
      }

      /**
       * Calculate the sensitivity RHS for the requested variables
       * (initials and/or parameters). Function signature is
       * pre-defined by CVODES library.
       *
       * @param[in] M number of parameters for which sensitivites are calculated
       * @param[in] t time
       * @param[in] y state of the base ODE system
       * @param[in] ydot state of the RHS of the ODE system
       * @param[in] yS array of M N_Vectors of size N, i.e. state of sensitivity
       * RHS
       * @param[out] ySdot array of M N_Vectors of size N of the sensitivity RHS
       */
      void rhs_sens(int M, realtype t,
                    double y[], double ydot[],
                    N_Vector *yS, N_Vector *ySdot) const {
        const std::vector<double> y_vec(y, y + N_);

        Eigen::VectorXd fy(N_);
        Eigen::MatrixXd Jy(N_, N_);

        if (param_var_) {
          Eigen::MatrixXd Jtheta(N_, M_);
          ode_model_.jacobian_SP(t, y_vec, fy, Jy, Jtheta);

          for (size_t m = 0; m < M_; m++) {
            // map NV_Vector to Eigen facilities
            Eigen::Map<Eigen::VectorXd>
              yS_eig(NV_DATA_S(yS[param_var_ind_ + m]), N_);
            Eigen::Map<Eigen::VectorXd>
              ySdot_eig(NV_DATA_S(ySdot[param_var_ind_ + m]), N_);

            ySdot_eig = Jy * yS_eig + Jtheta.col(m);
          }
        } else {
          ode_model_.jacobian_S(t, y_vec, fy, Jy);
        }

        if (initial_var_) {
          for (size_t m = 0; m < N_; m++) {
            // map NV_Vector to Eigen facilities
            Eigen::Map<Eigen::VectorXd> yS_eig(NV_DATA_S(yS[m]), N_);
            Eigen::Map<Eigen::VectorXd> ySdot_eig(NV_DATA_S(ySdot[m]), N_);

            ySdot_eig = Jy * yS_eig;
          }
        }
      }
      
      // Static wrapper for CVode callback
      static int ode_rhs_sens(int Ns, realtype t,
                              N_Vector y, N_Vector ydot,
                              N_Vector *yS, N_Vector *ySdot, void *user_data,
                              N_Vector tmp1, N_Vector tmp2) {
        ode const* explicit_ode = reinterpret_cast<ode const*>(user_data);
        explicit_ode->rhs_sens(Ns, t,
                               NV_DATA_S(y),  NV_DATA_S(ydot),
                               yS, ySdot);
        return 0;
      }

      /**
       * Integrate ODE and exract solution at time-points ts.
       *
       * @param[in] ts vector of time-points to evaluate
       * @param[out] y_coupled integrated function
       */
      void integrate_times(const std::vector<double>& ts,
                           std::vector<std::vector<double> >& y_coupled) const {
        double t_init = t0_;
        for (size_t n = 0; n < ts.size(); ++n) {
          double t_final = ts[n];
          if (t_final != t_init)
            check_flag(CVode(cvode_mem_, t_final, cvode_state_,
                             &t_init, CV_NORMAL),
                       "CVode");
          std::copy(state_.begin(), state_.end(), y_coupled[n].begin());
          if (S_ > 0) {
            check_flag(CVodeGetSens(cvode_mem_, &t_init, cvode_state_sens_),
                       "CVodeGetSens");
            for (size_t s = 0; s < S_; s++) {
              std::copy(NV_DATA_S(cvode_state_sens_[s]),
                        NV_DATA_S(cvode_state_sens_[s]) + N_,
                        y_coupled[n].begin() + N_ + s * N_);
            }
          }
          t_init = t_final;
        }
      }

      // Static wrapper for CVode callback
      static int dense_jacobian(long int N,  // NOLINT(runtime/int)
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *J_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode const* explicit_ode = reinterpret_cast<ode const*>(J_data);
        return explicit_ode->dense_jacobian(NV_DATA_S(y), J, t);
      }

      /**
       * Calculate the Jacobian of the ODE RHS wrt to states
       *
       * @param[in] y state of the base ODE system
       * @param[out] J Jacobian wrt to states y of ODE RHS
       * @param[in] t time
       */
      int dense_jacobian(const double* y, DlsMat J, double t) const {
        const std::vector<double> y_vec(y, y + N_);

        Eigen::VectorXd fy(N_);
        // Eigen and CVODES use column major addressing
        Eigen::Map<Eigen::MatrixXd> Jy_map(J->data, N_, N_);

        ode_model_.jacobian_S(t, y_vec, fy, Jy_map);

        return 0;
      }

      /**
       * Return the size of system (length of output vector at each
       * time-point).
       *
       * @return size of the theta vector
       */
      int size() const {
        return size_;
      }

      /**
       * Return the initial state of the coupled system. Note: Could
       * be removed, only here because tests expect it.
       *
       * @return initial state of coupled system
       */
      std::vector<double> initial_state() const {
        std::vector<double> state(size_, 0.0);
        std::copy(y0_dbl_.begin(), y0_dbl_.end(), state.begin());
        if (initial_var_) {
          for (size_t n = 0; n < N_; n++)
            state[N_ + n * N_ + n] = 1.0;
        }
        return state;
      }
    };

  }  // math
}  // stan
#endif
