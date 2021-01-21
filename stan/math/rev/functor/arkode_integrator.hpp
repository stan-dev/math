#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_ARKODE_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_ARKODE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/coupled_ode_system.hpp>
#include <stan/math/rev/functor/ode_store_sensitivities.hpp>
#include <stan/math/rev/functor/sundials_check.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <arkode/arkode.h>
#include <arkode/arkode_erkstep.h>
#include <arkode/arkode_butcher_erk.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Integrator interface for ARKode' ODE solvers
 *
 * @tparam scheme_id Butcher ID of ERK schemes.
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_param Type of scalars for parameters
 * @tparam T_t0 Type of scalar of initial time point
 * @tparam T_ts Type of time-points where ODE solution is returned
 */
template <int scheme_id, typename F, typename T_y0, typename T_t0,
          typename T_ts, typename... T_Args>
class arkode_integrator {
  using T_Return = return_type_t<T_y0, T_t0, T_ts, T_Args...>;
  using T_y0_t0 = return_type_t<T_y0, T_t0>;

  const F& f_;
  const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> y0_;
  const T_t0 t0_;
  const std::vector<T_ts>& ts_;
  std::tuple<const T_Args&...> args_tuple_;
  std::tuple<plain_type_t<decltype(value_of(std::declval<const T_Args&>()))>...>
      value_of_args_tuple_;
  const size_t N_;
  std::ostream* msgs_;
  double rtol_;
  double atol_;
  long int max_num_steps_;

  const size_t num_y0_vars_;
  const size_t num_args_vars_;

  coupled_ode_system<F, T_y0_t0, T_Args...> coupled_ode_;

  std::vector<double> coupled_state_;
  N_Vector nv_y_;

  void* mem_;

  /**
   * Implements the function of type CVRhsFn which is the user-defined
   * ODE RHS passed to ARKode.
   */
  static int arkode_rhs(realtype t, N_Vector y, N_Vector ydot,
                        void* user_data) {
    arkode_integrator* integrator = static_cast<arkode_integrator*>(user_data);
    integrator->rhs(t, y, ydot);
    return 0;
  }

  /**
   * Calculates the RHS of the sensitivity ODE system which
   * corresponds to the coupled ode system from which the first N
   * states are omitted, since the first N states are the ODE RHS
   * which ARKode separates from the main ODE RHS.
   */
  inline void rhs(double t, N_Vector y, N_Vector ydot) {
    size_t n = coupled_state_.size();
    std::vector<double> z(n);
    std::vector<double> dz_dt(n);
    std::copy(NV_DATA_S(y), NV_DATA_S(y) + n, z.begin());
    coupled_ode_(z, dz_dt, t);
    std::copy(dz_dt.begin(), dz_dt.end(), NV_DATA_S(ydot));
  }

 public:
  /**
   * Construct cvodes_integrator object
   *
   * @param function_name Calling function name (for printing debugging
   * messages)
   * @param f Right hand side of the ODE
   * @param y0 Initial state
   * @param t0 Initial time
   * @param ts Times at which to solve the ODE at. All values must be sorted and
   *   not less than t0.
   * @param relative_tolerance Relative tolerance passed to ARKode
   * @param absolute_tolerance Absolute tolerance passed to ARKode
   * @param max_num_steps Upper limit on the number of integration steps to
   *   take between each output (error if exceeded)
   * @param[in, out] msgs the print stream for warning messages
   * @param args Extra arguments passed unmodified through to ODE right hand
   * side function
   * @return Solution to ODE at times \p ts
   * @return a vector of states, each state being a vector of the
   *   same size as the state variable, corresponding to a time in ts.
   * @throw <code>std::domain_error</code> if y0, t0, ts, theta, x are not
   *   finite, all elements of ts are not greater than t0, or ts is not
   *   sorted in strictly increasing order.
   * @throw <code>std::invalid_argument</code> if arguments are the wrong
   *   size or tolerances or max_num_steps are out of range.
   */
  template <require_eigen_col_vector_t<T_y0>* = nullptr>
  arkode_integrator(const char* function_name, const F& f, const T_y0& y0,
                    const T_t0& t0, const std::vector<T_ts>& ts,
                    double relative_tolerance, double absolute_tolerance,
                    long int max_num_steps, std::ostream* msgs,
                    const T_Args&... args)
      : f_(f),
        y0_(y0.template cast<T_y0_t0>()),
        t0_(t0),
        ts_(ts),
        args_tuple_(args...),
        value_of_args_tuple_(value_of(args)...),
        N_(y0.size()),
        msgs_(msgs),
        rtol_(relative_tolerance),
        atol_(absolute_tolerance),
        max_num_steps_(max_num_steps),
        num_y0_vars_(count_vars(y0_)),
        num_args_vars_(count_vars(args...)),
        coupled_ode_(f, y0_, msgs, args...),
        coupled_state_(coupled_ode_.initial_state()),
        nv_y_(N_VMake_Serial(coupled_ode_.size(), &coupled_state_[0])),
        mem_(ERKStepCreate(arkode_rhs, value_of(t0_), nv_y_)) {
    check_finite(function_name, "initial state", y0_);
    check_finite(function_name, "initial time", t0_);
    check_finite(function_name, "times", ts_);

    // Code from: https://stackoverflow.com/a/17340003 . Should probably do
    // something better
    apply(
        [&](auto&&... args) {
          std::vector<int> unused_temp{
              0, (check_finite(function_name, "ode parameters and data", args),
                  0)...};
        },
        args_tuple_);

    check_nonzero_size(function_name, "times", ts_);
    check_nonzero_size(function_name, "initial state", y0_);
    check_sorted(function_name, "times", ts_);
    check_less(function_name, "initial time", t0_, ts_[0]);
    check_positive_finite(function_name, "relative_tolerance", rtol_);
    check_positive_finite(function_name, "absolute_tolerance", atol_);
    check_positive(function_name, "max_num_steps", max_num_steps_);

    if (mem_ == nullptr) {
      throw std::runtime_error("ERKStepCreate failed to allocate memory");
    }
  }

  ~arkode_integrator() {
    ERKStepFree(&mem_);
    N_VDestroy_Serial(nv_y_);
  }

  /**
   * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
   * times, { t1, t2, t3, ... } using the selected RK solver in ARKode.
   *
   * @return std::vector of Eigen::Matrix of the states of the ODE, one for each
   *   solution time (excluding the initial state)
   */
  std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> operator()() {
    std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> y;
    double t_init = value_of(t0_);

    CHECK_SUNDIALS_CALL(ERKStepSStolerances(mem_, rtol_, atol_));
    CHECK_SUNDIALS_CALL(ERKStepSetMaxNumSteps(mem_, max_num_steps_));
    CHECK_SUNDIALS_CALL(ERKStepSetTableNum(mem_, scheme_id));
    CHECK_SUNDIALS_CALL(ERKStepSetAdaptivityMethod(mem_, 2, 1, 0, 0));
    // CHECK_SUNDIALS_CALL(ERKStepSetInitStep(mem_, 0.1));
    ERKStepSetMaxEFailGrowth(mem_, 0.50);
    ERKStepSetMinReduction(mem_, 0.30);
    ERKStepSetFixedStepBounds(mem_, 1.0, 1.1);
    CHECK_SUNDIALS_CALL(ERKStepSetSafetyFactor(mem_, 0.9));
    CHECK_SUNDIALS_CALL(ERKStepSetUserData(mem_, reinterpret_cast<void*>(this)));

    for (size_t n = 0; n < ts_.size(); ++n) {
      double t_final = value_of(ts_[n]);

      if (t_final != t_init) {
        CHECK_SUNDIALS_CALL(
            ERKStepEvolve(mem_, t_final, nv_y_, &t_init, ARK_NORMAL));
      }

      y.emplace_back(apply(
          [&](auto&&... args) {
            return ode_store_sensitivities(f_, coupled_state_, y0_, t0_, ts_[n],
                                           msgs_, args...);
          },
          args_tuple_));

      // t_init = t_final;
    }

  /* Print some final statistics */
    static long int nst = 0, nst_a =0, nfe=0, netf=0, logging = false;
    int flag;
    // flag = ERKStepGetNumSteps(mem_, &nst);
    // flag = ERKStepGetNumStepAttempts(mem_, &nst_a);
    // flag = ERKStepGetNumRhsEvals(mem_, &nfe);
    // flag = ERKStepGetNumErrTestFails(mem_, &netf);

    // printf("\nFinal Solver Statistics:\n");
    // printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
    // printf("   Total RHS evals = %li\n", nfe);
    // printf("   Total number of error test failures = %li\n\n", netf);

    return y;
  }
};  // arkode integrator

}  // namespace math
}  // namespace stan
#endif
