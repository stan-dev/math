#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_CVODE_HPP

#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_band.h>

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
    private:
      void* cvode_mem_;
      vector<double> state_;
      N_Vector cvode_state_;

    public:
      typedef coupled_ode_system_cvode<F, T1, T2> ode;

      coupled_ode_system_cvode::~coupled_ode_system_cvode() {
        N_VDestroy_Serial(cvode_state_);
        CVodeFree(&cvode_mem_);
      }

      void check_flag(int flag, std::string func_name) {
        std::ostringstream ss;
        ss << func_name << " failed with error flag " << flag;
        throw std::runtime_error(ss.str());
      }

      // Establish all CVode preliminaries that
      // do not depend on the system size
      void instantiate_cvode_mem(double rel_tol,
                                 double abs_tol,
                                 long int max_num_steps,
                                 double init_step = 0
                                 long int max_err_test_fails = 20,
                                 long int max_conv_fails = 50) {
        // Instantiate CVode memory
        cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);
        if (cvode_mem == 0)
          throw std::runtime_error("CVodeCreate failed to allocate memory");

        // Forward CVode errors to noop error handler
        CVodeSetErrHandlerFn(cvode_mem_, silent_err_handler, 0);

        // Assign pointer to this as user data
        check_flag(CVodeSetUserData(cvode_mem, (void*)(this)),
                   "CVodeSetUserData");

        // Initialize solver parameters
        check_flag(CVodeSStolerances(cvode_mem, rel_tol, abs_tol),
                   "CVodeSStolerances");

        check_flag(CVodeSetMaxNumSteps(cvode_mem, max_num_steps),
                   "CVodeSetMaxNumSteps");

        check_flag(CVodeSetInitStep(cvode_mem, init_step),
                   "CVodeSetInitStep");

        check_flag(CVodeSetMaxErrTestFails(cvode_mem, max_err_test_fails),
                   "CVodeSetMaxErrTestFails");

        check_flag(CVodeSetMaxConvFails(cvode_mem, max_conv_fails),
                   "CVodeSetMaxConvFails");
      }

      // Forward the CVode array-based call to the ODE RHS
      // to the Stan vector-based call
      void operator()(const double y[], double dy_dy[], double t) {
        std::vector<double> y_(&y[0], &y[N_]);
        std::vector<double> dy_dt_;
        (*this)(y_, dy_dt_, t);
        std::copy(dy_dt_.begin(), dy_dt_.end(), &dy_dy[0]);
      }

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
                           std::vector<std::vector<double> >
                             y_coupled(ts_vec.size())) {
        check_flag(CVodeInit(cvode_mem_, &ode::ode_rhs, ts[0], cvode_state_),
                   "CVodeInit");
        y_coupled[0] = state_;
        for (int n = 1; n < ts.size(); ++n) {
          if (ts[n] != ts[n - 1])
            check_flag(CVode(cvode_mem_, ts[n], cvode_state_,
                             &(ts[n - 1]), CV_NORMAL),
                       "CVode");
          y_coupled[n] = state_;
        }
      }
    };

    /**
     * The coupled ode system for known initial values and known
     * parameters.  No sensitivities are needed and CVode uses
     * the full Jacobian to solve the system.
     *
     * @tparam F type of system function for the base ODE system.
     */
    template <typename F>
    class coupled_ode_system_cvode<F, double, double>:
      public coupled_ode_system<F, double, double> {
    public:
      // Static wrapper for CVode callback
      int dense_jacobian(long int N,
                         realtype t, N_Vector y, N_Vector fy,
                         DlsMat J, void *J_data,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
        ode* explicit_ode = static_cast<ode*>(J_data);
        return explicit_ode->dense_jacobian(NV_DATA_S(y), J->cols, t);
      }
      /**
       * Construct the coupled ODE system from the base system
       * function, initial state, parameters, data and a stream for
       * messages.
       *
       * @param[in] f base ode system functor.
       * @param[in] y0 initial state of the base ode.
       * @param[in] theta parameters of the base ode.
       * @param[in] x  real data.
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
        : coupled_ode_system(f, y0, theta, x, x_int, msgs),
          cvode_mem_(NULL),
          state_(size()),
          cvode_state_(size(), &state_[0]) {
        instantiate_cvode_mem(rel_tol, abs_tol, max_num_steps);
        check_flag(CVDense(cvode_mem_, size()), "CVDense");
        check_flag(CVDlsSetDenseJacFn(cvode_mem_, &ode::dense_jacobian),
                   "CVDlsSetDenseJacFn");
      }

      /**
       * Calculates the derivative of the coupled ode system with
       * respect to the specified state at the specified time using
       * the system state function.
       *
       * The derivative vector created is the same length as the
       * length as the state vector.
       *
       * @param[in] y current state of the coupled ode.
       * @param[out] dy_dt populated with derivatives of the coupled
       * system evaluated at specified state and time.
       * @param[in] t time.
       * @throw exception if the system function does not return
       * a derivative vector of the same size as the state vector.
       */
      void operator()(const double y[], double dy_dy[], double t) {
        std::vector<double> y_(&y[0], &y[size()]);
        std::vector<double> dy_dt_;
        (*this)(y_, dy_dt_, t);
        std::copy(dy_dt_.begin(), dy_dt_.end(), &dy_dy[0]);
      }

      // J[j][i] = d(ydot[i])/d(y[j])
      // J always initialized to zero by CVode
      void dense_jacobian(const double y[], double *J[], double t) {
        // compute Jacobian here!
      }

    };

  }  // math
}  // stan

/*
// J_data is a pointer to the ODE functor
extern "C"
int banded_jacobian(long int N, long int m_upper, long int m_lower,
                    realtype t, N_Vector y, N_Vector fy,
                    DlsMat J, void *J_data,
                    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  ode* explicit_ode = (ode*)J_data;
  return explicit_ode->banded_jacobian(t, NV_DATA_S(y), J->cols, J->s_mu);
}

int banded_jacobian(const double y[], double* J[], long int s_mu,
                    double t)

                    flag = CVBand(cvode_mem, n_total_state, n_state - 1, n_state - 1);
                    if (flag < 0) {
                      std::ostringstream ss;
                      ss << "CVBand failed with error flag " << flag;
                      throw std::runtime_error(ss.str());
                    }

                    flag = CVDlsSetBandJacFn(cvode_mem, banded_jacobian);
                    if (flag < 0) {
                      std::ostringstream ss;
                      ss << "CVDlsSetBandJacFn failed with error flag " << flag;
                      throw std::runtime_error(ss.str());
                    }

*/
#endif
