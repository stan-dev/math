#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP

#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/rev/core.hpp>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>

namespace stan {
  namespace math {

    // Noop error handler to silence CVode error output
    extern "C"
    void silent_err_handler(int error_code, const char *module,
                            const char *function, char *msg, void *eh_data) {
    }

    /**
     * Base template class for a coupled ordinary differential equation
     * system, which adds sensitivities to the base system, with additional
     * structure supporting CVode.
     *
     * This template class declaration should not be instantiated
     * directly --- it is just here to serve as a base for its
     * specializations, some of which are defined in namespace
     * <code>stan::agrad</code>.
     *
     * @tparam F the functor for the base ode system
     * @tparam T1 type of the initial state
     * @tparam T2 type of the parameters
     */
    template <typename F, typename T1, typename T2>
    class coupled_ode_system_cvode: public coupled_ode_system<F, T1, T2> {
    };

    /**
     * The coupled ODE system for known initial values and known
     * parameters.  No sensitivities are needed and CVode uses
     * the full Jacobian to solve the potentially-stiff system.
     *
     * @tparam F type of system function for the base ODE system.
     */
    template <typename F>
    class coupled_ode_system_cvode<F, double, double>:
      public coupled_ode_system<F, double, double> {
    private:
      void* cvode_mem_;
      std::vector<double> state_;
      N_Vector cvode_state_;

      typedef coupled_ode_system_cvode<F, double, double> ode;
      typedef coupled_ode_system<F, double, double> base_ode;

      void check_flag_(int flag, std::string func_name) {
        if (flag < 0) {
          std::ostringstream ss;
          ss << func_name << " failed with error flag " << flag;
          throw std::runtime_error(ss.str());
        }
      }

    public:
      /**
       * Construct the coupled ODE system from the base system
       * function, initial state, parameters, data and a stream for
       * messages.
       *
       * @param[in] f base ode system functor.
       * @param[in] y0 initial state of the base ode.
       * @param[in] theta parameters of the base ode.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in] max_num_steps Maximum number of solver steps.
       * @param[in] rel_tol Relative tolerance of solver.
       * @param[in] abs_tol Absolute tolerance of solver.
       * @param[in, out] msgs print stream.
       */
      coupled_ode_system_cvode(const F& f,
                               const std::vector<double>& y0,
                               const std::vector<double>& theta,
                               const std::vector<double>& x,
                               const std::vector<int>& x_int,
                               double rel_tol,
                               double abs_tol,
                               long int max_num_steps,
                               std::ostream* msgs)
        : coupled_ode_system<F, double, double>(f, y0, theta, x, x_int, msgs),
          cvode_mem_(NULL),
          state_(this->N_),
          cvode_state_(N_VMake_Serial(this->N_, &state_[0])) {
        // Instantiate CVode memory
        cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);
        if (cvode_mem_ == 0)
          throw std::runtime_error("CVodeCreate failed to allocate memory");

        // Forward CVode errors to noop error handler
        CVodeSetErrHandlerFn(cvode_mem_, silent_err_handler, 0);

        // Assign pointer to this as user data
        check_flag_(CVodeSetUserData(cvode_mem_, (void*)(this)),
                    "CVodeSetUserData");

        // Initialize solver parameters
        check_flag_(CVodeSStolerances(cvode_mem_, rel_tol, abs_tol),
                    "CVodeSStolerances");

        check_flag_(CVodeSetMaxNumSteps(cvode_mem_, max_num_steps),
                    "CVodeSetMaxNumSteps");

        double init_step = 0;
        check_flag_(CVodeSetInitStep(cvode_mem_, init_step),
                    "CVodeSetInitStep");

        long int max_err_test_fails = 20;
        check_flag_(CVodeSetMaxErrTestFails(cvode_mem_, max_err_test_fails),
                    "CVodeSetMaxErrTestFails");

        long int max_conv_fails = 50;
        check_flag_(CVodeSetMaxConvFails(cvode_mem_, max_conv_fails),
                    "CVodeSetMaxConvFails");

        check_flag_(CVDense(cvode_mem_, this->N_), "CVDense");
        check_flag_(CVDlsSetDenseJacFn(cvode_mem_, &ode::dense_jacobian),
               "CVDlsSetDenseJacFn");
      }

      ~coupled_ode_system_cvode() {
        N_VDestroy_Serial(cvode_state_);
        CVodeFree(&cvode_mem_);
      }

      using coupled_ode_system<F, double, double>::operator();

      // Static wrapper for CVode callback
      static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* f_data) {
        static_cast<ode*>(f_data)->operator()(NV_DATA_S(y), NV_DATA_S(ydot), t);
      }
      /*
      static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* f_data) {
        coupled_ode_system_cvode<F, T1, T2>* explicit_ode
          = static_cast<coupled_ode_system_cvode<F, T1, T2>*>(f_data);
        return explicit_ode->operator()(NV_DATA_S(y), NV_DATA_S(ydot), t);
      }*/

      void integrate_times(const std::vector<double>& ts,
                           std::vector<std::vector<double> > y_coupled) {
        check_flag_(CVodeInit(cvode_mem_, &ode::ode_rhs, ts[0], cvode_state_),
                    "CVodeInit");
        y_coupled[0] = state_; // <|--- what is the optimal way to copy vectors?
        for (int n = 1; n < ts.size(); ++n) {
          if (ts[n] != ts[n - 1])
            check_flag_(CVode(cvode_mem_, ts[n], cvode_state_,
                              &(ts[n - 1]), CV_NORMAL),
                        "CVode");
          y_coupled[n] = state_; // <|--- what is the optimal way to copy vectors?
        }
      }

      // Forward the CVode array-based call from
      // the ODE RHS to the Stan vector-based call
      void operator()(const double y[], double dy_dt[], double t) {
        std::vector<double> y_vec(this->N_);
        std::copy(y, y + this->N_, y_vec);

        std::vector<double> dy_dt_vec(this->N_);
        (*this)(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int dense_jacobian(long int N,
                                realtype t, N_Vector y, N_Vector fy,
                                DlsMat J, void *J_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode* explicit_ode = static_cast<ode*>(J_data);
        return explicit_ode->dense_jacobian(NV_DATA_S(y), J->cols, t);
      }

      // J[j][i] = d(ydot[i])/d(y[j])
      int dense_jacobian(const double y[], double *J[], double t) {
        using stan::math::var;

        std::vector<double> y_vec(this->N_);
        for (int n = 0; n < this->N_; ++n)
          y_vec[n] = y[n];
        //std::copy(y, y + this->N_, y_vec);

        std::vector<double> grad(this->N_);

        try {
          stan::math::start_nested();

          std::vector<var> y_vars(y_vec.begin(), y_vec.end());
          std::vector<var> dy_dt_vars = this->f_(t, y_vars, this->theta_dbl_,
                                                 this->x_, this->x_int_,
                                                 this->msgs_);

          for (size_t i = 0; i < this->N_; i++) {
            set_zero_all_adjoints_nested();
            dy_dt_vars[i].grad(y_vars, grad);
            for (size_t j = 0; j < this->N_; ++j)
              J[j][i] = grad[j];
          }
        } catch (const std::exception& e) {
          stan::math::recover_memory_nested();
          throw;
        }
        stan::math::recover_memory_nested();
        return 0;
      }
    };
  }  // math
}  // stan

#endif
