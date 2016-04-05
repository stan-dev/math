#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_CVODES_INTEGRATOR_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_CVODES_INTEGRATOR_HPP

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/arr/functor/ode_model.hpp>

#include <vector>
#include <algorithm>
#include <string>
#include <boost/type_traits/is_same.hpp>

namespace stan {
  namespace math {

    // Noop error handler to silence CVodes error output
    extern "C"
    void silent_err_handler(int error_code, const char *module,
                            const char *function, char *msg, void *eh_data) {
    }

    template <typename F, typename T1, typename T2>
    class cvodes_integrator {
    private:
      double t0_;
      const std::vector<double>& y0_dbl_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const size_t theta_var_ind_;
      void* cvode_mem_;
      std::vector<double> state_;
      N_Vector  cvode_state_;
      N_Vector *cvode_state_sens_;

      const ode_model<F> ode_model_;
      const size_t solver_;

      typedef cvodes_integrator<F, T1, T2> ode;

      void check_flag_(int flag, std::string func_name) {
        if (flag < 0) {
          std::ostringstream ss;
          ss << func_name << " failed with error flag " << flag;
          throw std::runtime_error(ss.str());
        }
      }

    public:
      // returned solution type
      typedef std::vector<typename stan::return_type<T1,T2>::type> state_t;
      typedef boost::is_same<T1, stan::math::var> initial_var;
      typedef boost::is_same<T2, stan::math::var> theta_var;

      /**
       * Construct the coupled ODE system from the base system
       * function, initial state, parameters, data and a stream for
       * messages.
       *
       * @param[in] f base ode system functor.
       * @param[in] y0 initial state of the base ode.
       * @param[in] t0 initial time of the base ode.
       * @param[in] theta parameters of the base ode.
       * @param[in] x real data.
       * @param[in] x_int integer data.
       * @param[in] rel_tol Relative tolerance of solver.
       * @param[in] abs_tol Absolute tolerance of solver.
       * @param[in] max_num_steps Maximum number of solver steps.
       * @param[in] solver used solver (0=non-stiff, 1=stiff, 2=stiff
       * with STALD)
       * @param[in, out] msgs print stream.
       */
      cvodes_integrator(const F& f,
			const std::vector<double>& y0,
			double t0,
			const std::vector<double>& theta,
			const std::vector<double>& x,
			const std::vector<int>& x_int,
			double rel_tol,
			double abs_tol,
			long int max_num_steps,  // NOLINT(runtime/int)
			size_t solver,
			std::ostream* msgs)
	: t0_(t0),
	  y0_dbl_(y0),
	  N_(y0.size()),
	  M_(theta.size()),
	  S_((initial_var::value ? N_ : 0) + (theta_var::value ? M_ : 0)),
	  size_(N_ * (S_+1)),
	  theta_var_ind_(initial_var::value ? N_ : 0),
	  cvode_mem_(NULL),
	  state_(y0),
	  cvode_state_(N_VMake_Serial(this->N_,
				      &state_[0])),
	  ode_model_(f, theta, x, x_int, msgs),
	  solver_(solver)
      {
        // Instantiate CVode memory
	if(solver_ == 0)
	  cvode_mem_ = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	if(solver_ == 1 || solver_ == 2)
	  cvode_mem_ = CVodeCreate(CV_BDF, CV_NEWTON);

        if (cvode_mem_ == 0)
          throw std::runtime_error("CVodeCreate failed to allocate memory");

        // Forward CVode errors to noop error handler
        CVodeSetErrHandlerFn(cvode_mem_, silent_err_handler, 0);

        check_flag_(CVodeInit(cvode_mem_, &ode::ode_rhs, t0_, cvode_state_),
                    "CVodeInit");

        // Assign pointer to this as user data
        check_flag_(CVodeSetUserData(cvode_mem_,
                                     reinterpret_cast<void*>(this)),
                    "CVodeSetUserData");

        // Initialize solver parameters
        check_flag_(CVodeSStolerances(cvode_mem_, rel_tol, abs_tol),
                    "CVodeSStolerances");

        check_flag_(CVodeSetMaxNumSteps(cvode_mem_, max_num_steps),
                    "CVodeSetMaxNumSteps");

        double init_step = 0;
        check_flag_(CVodeSetInitStep(cvode_mem_, init_step),
                    "CVodeSetInitStep");

        long int max_err_test_fails = 20;  // NOLINT(runtime/int)
        check_flag_(CVodeSetMaxErrTestFails(cvode_mem_, max_err_test_fails),
                    "CVodeSetMaxErrTestFails");

        long int max_conv_fails = 50;  // NOLINT(runtime/int)
        check_flag_(CVodeSetMaxConvFails(cvode_mem_, max_conv_fails),
                    "CVodeSetMaxConvFails");

	// enable stability limit detection for solver_==2
	if(solver_ == 2)
          check_flag_(CVodeSetStabLimDet(cvode_mem_, 1), "CVodeSetstabLimDet");

	if(solver_ == 1 || solver_ == 2) {
	  // for the stiff solver we need to reserve memory and
	  // provide a Jacobian
	  check_flag_(CVDense(cvode_mem_, this->N_), "CVDense");
	  check_flag_(CVDlsSetDenseJacFn(cvode_mem_, &ode::dense_jacobian),
		      "CVDlsSetDenseJacFn");
	}

	// initialize sensitivity system for parameters as needed
	if(S_ > 0) {
	  cvode_state_sens_ = N_VCloneVectorArray_Serial(S_, cvode_state_);
	  for(size_t s=0; s < S_ ; s++) N_VConst(RCONST(0.0), cvode_state_sens_[s]);
	  // for the case with varying initials, the first N_
	  // sensitivity systems correspond to the initials which have
	  // as initial the identity matrix
	  if(initial_var::value) {
	    for(size_t n=0; n < N_; n++) {
	      NV_Ith_S(cvode_state_sens_[n],n) = 1.0;
	    }
	  }

	  check_flag_(CVodeSensInit(cvode_mem_, (int)S_, CV_STAGGERED, &ode::ode_rhs_sens, cvode_state_sens_),
		      "CVodeSensInit");

	  check_flag_(CVodeSensEEtolerances(cvode_mem_),
		      "CVodeSensEEtolerances");
	}

      }

      ~cvodes_integrator() {
        // N_VDestroy_Serial should noop because
        // N_VMake_Serial sets own_data to false
        N_VDestroy_Serial(cvode_state_);
	if(S_ > 0)
	  N_VDestroyVectorArray_Serial(cvode_state_sens_, S_);
        CVodeFree(&cvode_mem_);
      }

      // Forward the CVode array-based call from
      // the ODE RHS to the Stan vector-based call
      void rhs(const double* y, double* dy_dt, double t) {
        const std::vector<double> y_vec(y, y + this->N_);

        std::vector<double> dy_dt_vec(this->N_);
	ode_model_(y_vec, dy_dt_vec, t);

        std::copy(dy_dt_vec.begin(), dy_dt_vec.end(), dy_dt);
      }

      // Static wrapper for CVode callback
      static int ode_rhs(double t, N_Vector y, N_Vector ydot, void* f_data) {
        ode* explicit_ode = reinterpret_cast<ode*>(f_data);
        explicit_ode->rhs(NV_DATA_S(y), NV_DATA_S(ydot), t);
        return 0;
      }

      void rhs_sens(int M, realtype t,
		    double y[], double ydot[],
		    N_Vector *yS, N_Vector *ySdot) {
        const std::vector<double> y_vec(y, y + N_);

	Eigen::VectorXd fy(N_);
	Eigen::MatrixXd Jy(N_,N_);
	Eigen::MatrixXd Jtheta(N_,M_);

	// do AD
	if(theta_var::value)
	  ode_model_.jacobian_SP(t, y_vec, fy, Jy, Jtheta);
	else
	  ode_model_.jacobian_S(t, y_vec, fy, Jy);

	if(initial_var::value) {
	  for(size_t m = 0; m < N_; m++) {
	    // map NV_Vector to Eigen facilities
	    Eigen::Map<Eigen::VectorXd>    yS_eig(NV_DATA_S(   yS[m]), N_);
	    Eigen::Map<Eigen::VectorXd> ySdot_eig(NV_DATA_S(ySdot[m]), N_);

	    ySdot_eig = Jy * yS_eig;
	  }
	}

	if(theta_var::value) {
	  for(size_t m = 0; m < M_; m++) {
	    // map NV_Vector to Eigen facilities
	    Eigen::Map<Eigen::VectorXd>    yS_eig(NV_DATA_S(   yS[theta_var_ind_ + m]), N_);
	    Eigen::Map<Eigen::VectorXd> ySdot_eig(NV_DATA_S(ySdot[theta_var_ind_ + m]), N_);

	    ySdot_eig = Jy * yS_eig + Jtheta.col(m);
	  }
	}

      }

      static int ode_rhs_sens(int Ns, realtype t,
			      N_Vector y, N_Vector ydot,
			      N_Vector *yS, N_Vector *ySdot, void *user_data,
			      N_Vector tmp1, N_Vector tmp2) {
        ode* explicit_ode = reinterpret_cast<ode*>(user_data);
        explicit_ode->rhs_sens(Ns, t,
			       NV_DATA_S(y),  NV_DATA_S(ydot),
			       yS, ySdot);
        return 0;
      }

      void report_stats(long int& n_rhs, long int& n_rhs_sens) {
	check_flag_(CVodeGetNumRhsEvals(cvode_mem_, &n_rhs), "CVodeGetNumRhsEvals");
	n_rhs_sens = 0;
	if(S_ > 0) {
	  check_flag_(CVodeGetSensNumRhsEvals(cvode_mem_, &n_rhs_sens), "CVodeGetSensNumRhsEvals");
	}
      }

      void integrate_times(const std::vector<double>& ts,
                           std::vector<std::vector<double> >& y_coupled) {
        double t_init = t0_;
        for (size_t n = 0; n < ts.size(); ++n) {
          double t_final = ts[n];
          if (t_final != t_init)
            check_flag_(CVode(cvode_mem_, t_final, cvode_state_,
                              &t_init, CV_NORMAL),
                        "CVode");
          std::copy(state_.begin(), state_.end(), y_coupled[n].begin());
	  if(S_ > 0) {
	    check_flag_(CVodeGetSens(cvode_mem_, &t_init, cvode_state_sens_), "CVodeGetSens");
	    for(size_t s = 0; s < S_; s++) {
	      for(size_t i = 0; i < N_; i++) {
		y_coupled[n][N_ + s * N_ + i] = NV_Ith_S(cvode_state_sens_[s],i);
	      }
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
        ode* explicit_ode = reinterpret_cast<ode*>(J_data);
	return explicit_ode->dense_jacobian(NV_DATA_S(y), J, t);
      }

      int dense_jacobian(const double* y, DlsMat J, double t) const {
        const std::vector<double> y_vec(y, y + this->N_);

	Eigen::VectorXd fy(this->N_);
	// map Eigen matrix to memory chunk where the data is to be
	// stored (Eigen and CVODES use column major addressing)
	Eigen::Map<Eigen::MatrixXd> Jy_map(J->data, this->N_, this->N_);

	ode_model_.jacobian_S(t, y_vec, fy, Jy_map);

        return 0;
      }

      int size() const {
        return size_;
      }

      std::vector<double> initial_state() {
        std::vector<double> state(size_, 0.0);
        for (size_t n = 0; n < N_; n++)
          state[n] = y0_dbl_[n];
	if(initial_var::value) {
	  for (size_t n = 0; n < N_; n++)
	    state[N_ + n * N_ + n] = 1;
	}
        return state;
      }

    };

  }  // math
}  // stan
#endif
