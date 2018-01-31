#ifndef STAN_MATH_REV_MAT_FUNCTOR_IDAS_INTEGRATOR_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_IDAS_INTEGRATOR_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/idas_system.hpp>
// #include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
#include <idas/idas.h>
#include <idas/idas_direct.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <ostream>
#include <vector>


enum IDAS_SENSITIVITY { forward, adjoint };

namespace stan {
namespace math {

/**
 * Integrator interface for CVODES' ODE solvers (Adams & BDF
 * methods).
 */
  // template<typename S>
  class idas_integrator {
    
    template<typename DAE>
    void init_sensitivity(DAE& dae);

    template<typename F, typename TYY, typename TYP, typename TPAR>
    void init_sensitivity(idas_adjoint_system<F, TYY, TYP, TPAR>& dae){
      // TODO
    }

    template<typename DAE>
    void solve(DAE& dae, double t0,
               const std::vector<double>& ts,
               std::vector<std::vector<double> > res_yy,
               std::vector<Eigen::MatrixXd> res_ys);

    template<typename F, typename TYY, typename TYP, typename TPAR>
    void solve(idas_adjoint_system<F, TYY, TYP, TPAR>& dae, double t0,
               const std::vector<double>& ts,
               std::vector<std::vector<double> > res_yy,
               std::vector<Eigen::MatrixXd> res_ys){
      // TODO
    }


  public:
    idas_integrator() {}
    
    /**
     * Return the solutions for the specified DAE
     * given the specified initial state,
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
     * @param[in, out] msgs the print stream for warning messages.
     * @param[in] relative_tolerance relative tolerance passed to CVODE.
     * @param[in] absolute_tolerance absolute tolerance passed to CVODE.
     * @param[in] max_num_steps maximal number of admissable steps
     * between time-points
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template<typename DAE>
    typename DAE::return_type
    integrate(DAE& dae,
              double t0,
              const std::vector<double>& ts,
              double rtol,
              double atol,
              long int max_num_steps) {

      using Eigen::Matrix;
      using Eigen::MatrixXd;
      using Eigen::VectorXd;
      using Eigen::Dynamic;

      //   const auto caller {"idas_integrator"};
      //   check_finite(caller, "initial state", yy0);
      //   check_finite(caller, "derivative initial state", yp0);
      //   check_consistent_sizes(caller, "initial state", yy0, "derivative initial state", yp0);
      //   check_greater_or_equal(caller, "derivative-algebra id", dae_id, 0);
      //   check_less_or_equal(caller, "derivative-algebra id", dae_id, 1);
      //   check_finite(caller, "initial time", t0);
      //   check_finite(caller, "times", ts);
      //   check_finite(caller, "parameter vector", theta);
      //   check_finite(caller, "continuous data", x_r);
      //   check_nonzero_size(caller, "times", ts);
      //   check_less(caller, "initial time", t0, ts.front());
      // if (rtol <= 0) invalid_argument(caller, "relative_tolerance,",
      //                                 rtol, "", ", must be greater than 0");
      // if (atol <= 0) invalid_argument(caller, "absolute_tolerance,",
      //                                 atol, "", ", must be greater than 0");
      // if (max_num_steps <= 0)
      //   invalid_argument(caller, "max_num_steps,", max_num_steps,
      //                    "", ", must be greater than 0");

      auto mem = dae.mem();
      auto yy = dae.nv_yy();
      auto yp = dae.nv_yp();
      auto n = dae.n();

      std::vector<std::vector<double> > res_yy(ts.size(),
                                               std::vector<double>(dae.n_sys(), 0));
      std::vector<MatrixXd> res_ys(ts.size(),
                                   MatrixXd::Zero(n, dae.n_sens()));

      // void* user_data = static_cast<void*>(&dae);

      // auto callback=[](double t, N_Vector yy, N_Vector yp,
      //                      N_Vector rr, void *user_data) -> int {
      //   // return dae.residual(t, yy, yp, rr, user_data);
      //   DAE* dae = static_cast<DAE*>(user_data);

      // };

      // auto callback=[](void* arg){ // captureless
      //   (*static_cast<decltype(dae.residual)>(arg))();
      // };

      try {
        CHECK_IDAS_CALL(IDAInit(mem, dae.residual(), t0, yy, yp));
        CHECK_IDAS_CALL(IDADlsSetLinearSolver(mem, dae.linsol(), dae.jacobi()));
        CHECK_IDAS_CALL(IDASStolerances(mem, rtol, atol));
        CHECK_IDAS_CALL(IDASetMaxNumSteps(mem, max_num_steps));
        dae.cast_to_user_data();
        dae.set_consistent_ic(ts.front());
        init_sensitivity(dae);

        solve(dae, t0, ts, res_yy, res_ys);

      } catch (const std::exception& e) {
        throw;
      }

      // return dae.solution(res_yy, res_ys);

    }
  };  // idas integrator


  template<typename DAE>
  void idas_integrator::init_sensitivity(DAE& dae) {
    if (DAE::need_sens) {
      auto mem = dae.mem();
      auto yys = dae.nv_yys();
      auto yps = dae.nv_yps();
      auto n = dae.n();
      CHECK_IDAS_CALL(IDASensInit(mem, dae.n_sens(),
                                  IDA_SIMULTANEOUS,
                                  dae.sensitivity_residual(),
                                  yys,
                                  yps));
      CHECK_IDAS_CALL(IDASensEEtolerances(mem));
      dae.set_consistent_sens_ic();
      if (DAE::is_var_yy0) for (size_t i = 0; i < n; i++) NV_Ith_S(yys[i], i) = 1.0;
      if (DAE::is_var_yp0) for (size_t i = 0; i < n; i++) NV_Ith_S(yps[i+n], i) = 1.0;
    }
  }

  template<typename DAE>
  void idas_integrator::solve(DAE& dae, double t0,
                              const std::vector<double>& ts,
                              std::vector<std::vector<double> > res_yy,
                              std::vector<Eigen::MatrixXd> res_ys){
    double t1 = t0;
    size_t i = 0;
    auto mem = dae.mem();
    auto yy = dae.nv_yy();
    auto yp = dae.nv_yp();
    auto yys = dae.nv_yys();
    auto yps = dae.nv_yps();
    auto n = dae.n();

    std::for_each(ts.begin(), ts.end(), [&](double t2){
        CHECK_IDAS_CALL(IDASolve(mem, t2, &t1, yy, yp, IDA_NORMAL));
        std::copy(dae.yy().begin(), dae.yy().end(), res_yy[i].begin());
        if(DAE::need_sens) {
          CHECK_IDAS_CALL(IDAGetSens(mem, &t1, yys));
          for(int j=0; j<dae.n_sens(); ++j) {
            for(int k=0; k<n; ++k) {
              res_ys[i](k, j) = NV_Ith_S(yys[j], k);
            }
          }
        }
        ++i;
      });
  }

}  // namespace math
}  // namespace stan

#endif
