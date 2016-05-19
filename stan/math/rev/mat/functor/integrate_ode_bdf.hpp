#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_BDF_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_BDF_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/cvodes_utils.hpp>
#include <stan/math/rev/mat/functor/cvodes_ode_data.hpp>
#include <stan/math/rev/arr/fun/decouple_ode_states.hpp>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>

#include <ostream>
#include <vector>
#include <algorithm>

namespace stan {
  namespace math {
    /**
     * Return the solutions for the specified system of ordinary
     * differential equations given the specified initial state,
     * initial times, times of desired solution, and parameters and
     * data, writing error and warning messages to the specified
     * stream.
     *
     * This function is templated to allow the initial times to be
     * either data or autodiff variables and the parameters to be data
     * or autodiff variables.  The autodiff-based implementation for
     * reverse-mode are defined in namespace <code>stan::math</code>
     * and may be invoked via argument-dependent lookup by including
     * their headers.
     *
     * The solver used is based on the backward differentiation
     * formula which is an implicit numerical integration scheme
     * appropiate for stiff ODE systems.
     *
     * @tparam F type of ODE system function.
     * @tparam T_initial type of scalars for initial values.
     * @tparam T_param type of scalars for parameters.
     * @param[in] f functor for the base ordinary differential equation.
     * @param[in] y0 initial state.
     * @param[in] t0 initial time.
     * @param[in] ts times of the desired solutions, in strictly
     * increasing order, all greater than the initial time.
     * @param[in] theta parameter vector for the ODE.
     * @param[in] x continuous data vector for the ODE.
     * @param[in] x_int integer data vector for the ODE.
     * @param[in] rel_tol relative tolerance of solution
     * @param[in] abs_tol absolute tolerance of solution
     * @param[in] max_num_steps maximal number of admissable steps
     * between time-points
     * @param[in, out] msgs the print stream for warning messages.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename F, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    integrate_ode_bdf(const F& f,
                      const std::vector<T_initial>& y0,
                      const double t0,
                      const std::vector<double>& ts,
                      const std::vector<T_param>& theta,
                      const std::vector<double>& x,
                      const std::vector<int>& x_int,
                      double rel_tol = 1e-10,
                      double abs_tol = 1e-10,
                      long int max_num_steps = 1e8,  // NOLINT(runtime/int)
                      std::ostream* msgs = 0) {
      stan::math::check_finite("integrate_ode_bdf",
                               "initial state", y0);
      stan::math::check_finite("integrate_ode_bdf",
                               "initial time", t0);
      stan::math::check_finite("integrate_ode_bdf",
                               "times", ts);
      stan::math::check_finite("integrate_ode_bdf",
                               "parameter vector", theta);
      stan::math::check_finite("integrate_ode_bdf",
                               "continuous data", x);

      stan::math::check_nonzero_size("integrate_ode_bdf",
                                     "times", ts);
      stan::math::check_nonzero_size("integrate_ode_bdf",
                                     "initial state", y0);
      stan::math::check_ordered("integrate_ode_bdf",
                                "times", ts);
      stan::math::check_less("integrate_ode_bdf",
                             "initial time", t0, ts[0]);

      typedef stan::is_var<T_initial> initial_var;
      typedef stan::is_var<T_param> param_var;

      // number of states
      const size_t N = y0.size();
      // number of parameters
      const size_t M = theta.size();
      // number of sensitivities, can be 0+0, N+0, 0+M, N+M for
      // varying/constants initials + varying/constant parameters
      const size_t S = (initial_var::value ? N : 0)
        + (param_var::value ? M : 0);
      // total size of the coupled system, N*(S+1)
      const size_t size = N * (S + 1);
      std::vector<double> state(stan::math::value_of(y0));
      N_Vector  cvodes_state(N_VMake_Serial(N, &state[0]));
      N_Vector *cvodes_state_sens = NULL;

      typedef cvodes_ode_data<F, T_initial, T_param> ode_data;
      ode_data cvodes_data(f, y0, theta, x, x_int, msgs);

      void* cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);

      if (cvodes_mem == NULL)
        throw std::runtime_error("CVodeCreate failed to allocate memory");

      std::vector<std::vector<double> >
        y_coupled(ts.size(), std::vector<double>(size, 0));

      try {
        cvodes_check_flag(CVodeInit(cvodes_mem, &ode_data::ode_rhs,
                                    t0, cvodes_state),
                          "CVodeInit");

        // Assign pointer to this as user data
        cvodes_check_flag(CVodeSetUserData(cvodes_mem,
                                   reinterpret_cast<void*>(&cvodes_data)),
                          "CVodeSetUserData");

        cvodes_set_options(cvodes_mem, rel_tol, abs_tol, max_num_steps);

        // for the stiff solvers we need to reserve additional
        // memory and provide a Jacobian function call
        cvodes_check_flag(CVDense(cvodes_mem, N), "CVDense");
        cvodes_check_flag(CVDlsSetDenseJacFn(cvodes_mem,
                                             &ode_data::dense_jacobian),
                          "CVDlsSetDenseJacFn");

        // initialize forward sensitivity system of CVODES as needed
        if (S > 0) {
          cvodes_state_sens = N_VCloneVectorArray_Serial(S, cvodes_state);

          for (size_t s = 0; s < S; s++)
            N_VConst(RCONST(0.0), cvodes_state_sens[s]);

          // for the case with varying initials, the first N
          // sensitivity systems correspond to the initials which have
          // as initial the identity matrix
          if (initial_var::value) {
            for (size_t n = 0; n < N; n++) {
              NV_Ith_S(cvodes_state_sens[n], n) = 1.0;
            }
          }

          cvodes_check_flag(CVodeSensInit(cvodes_mem, static_cast<int>(S),
                                          CV_STAGGERED,
                                          &ode_data::ode_rhs_sens,
                                          cvodes_state_sens),
                            "CVodeSensInit");

          cvodes_check_flag(CVodeSensEEtolerances(cvodes_mem),
                            "CVodeSensEEtolerances");
        }

        double t_init = t0;
        for (size_t n = 0; n < ts.size(); ++n) {
          double t_final = ts[n];
          if (t_final != t_init)
            cvodes_check_flag(CVode(cvodes_mem, t_final, cvodes_state,
                                    &t_init, CV_NORMAL),
                              "CVode");
          std::copy(state.begin(), state.end(), y_coupled[n].begin());
          if (S > 0) {
            cvodes_check_flag(CVodeGetSens(cvodes_mem, &t_init,
                                           cvodes_state_sens),
                              "CVodeGetSens");
            for (size_t s = 0; s < S; s++) {
              std::copy(NV_DATA_S(cvodes_state_sens[s]),
                        NV_DATA_S(cvodes_state_sens[s]) + N,
                        y_coupled[n].begin() + N + s * N);
            }
          }
          t_init = t_final;
        }
      } catch (const std::exception& e) {
        N_VDestroy_Serial(cvodes_state);
        if (cvodes_state_sens != NULL)
          N_VDestroyVectorArray_Serial(cvodes_state_sens, S);
        CVodeFree(&cvodes_mem);
        throw;
      }

      N_VDestroy_Serial(cvodes_state);
      if (cvodes_state_sens != NULL)
        N_VDestroyVectorArray_Serial(cvodes_state_sens, S);
      CVodeFree(&cvodes_mem);

      return decouple_ode_states(y_coupled, y0, theta);
    }
  }  // math
}  // stan
#endif
