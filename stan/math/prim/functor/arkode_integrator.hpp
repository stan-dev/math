#ifndef STAN_MATH_PRIM_FUNCTOR_ARKODE_INTEGRATOR_HPP
#define STAN_MATH_PRIM_FUNCTOR_ARKODE_INTEGRATOR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <sundials/sundials_context.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Integrator interface for ARKODE's ERKStep explicit adaptive Runge-Kutta
 * solver, using one of ERKStep's built-in embedded Butcher tables.
 *
 * @tparam Table ID of the built-in ERKStep Butcher table to use, e.g.
 *   <code>ARKODE_CASH_KARP_6_4_5</code> or
 *   <code>ARKODE_DORMAND_PRINCE_7_4_5</code>
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_t0 Type of scalar of initial time point
 * @tparam T_ts Type of time-points where ODE solution is returned
 */
template <ARKODE_ERKTableID Table, typename F, typename T_y0, typename T_t0,
          typename T_ts, typename... Args>
class arkode_integrator {
  using T_Return = return_type_t<T_y0, T_t0, T_ts, Args...>;
  using T_y0_t0 = return_type_t<T_y0, T_t0>;

  const char* function_name_;
  sundials::Context sundials_context_;
  const F& f_;
  const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> y0_;
  const T_t0 t0_;
  const std::vector<T_ts>& ts_;
  std::tuple<const Args&...> args_tuple_;
  std::ostream* msgs_;
  double relative_tolerance_;
  double absolute_tolerance_;
  long int max_num_steps_;  // NOLINT(runtime/int)

  coupled_ode_system<F, T_y0_t0, Args...> coupled_system_;
  std::vector<double> coupled_state_;
  N_Vector nv_state_;

  static int erk_rhs(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    arkode_integrator* integrator = static_cast<arkode_integrator*>(user_data);
    integrator->rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return 0;
  }

  inline void rhs(double t, const double y[], double dy_dt[]) {
    std::vector<double> z(y, y + coupled_state_.size());
    std::vector<double> dz_dt;
    coupled_system_(z, dz_dt, t);
    std::copy(dz_dt.begin(), dz_dt.end(), dy_dt);
  }

 public:
  /**
   * Construct arkode_integrator object
   *
   * @param function_name Calling function name (for printing debugging
   * messages)
   * @param f Right hand side of the ODE
   * @param y0 Initial state
   * @param t0 Initial time
   * @param ts Times at which to solve the ODE at. All values must be sorted
   *   and greater than t0.
   * @param relative_tolerance Relative tolerance passed to ARKODE
   * @param absolute_tolerance Absolute tolerance passed to ARKODE
   * @param max_num_steps Upper limit on the number of integration steps to
   *   take between each output (error if exceeded)
   * @param[in, out] msgs the print stream for warning messages
   * @param args Extra arguments passed unmodified through to ODE right hand
   * side function
   */
  arkode_integrator(const char* function_name, const F& f, const T_y0& y0,
                    T_t0 t0, const std::vector<T_ts>& ts,
                    double relative_tolerance, double absolute_tolerance,
                    long int max_num_steps,  // NOLINT(runtime/int)
                    std::ostream* msgs, const Args&... args)
      : function_name_(function_name),
        sundials_context_(),
        f_(f),
        y0_(y0.template cast<T_y0_t0>()),
        t0_(t0),
        ts_(ts),
        args_tuple_(args...),
        msgs_(msgs),
        relative_tolerance_(relative_tolerance),
        absolute_tolerance_(absolute_tolerance),
        max_num_steps_(max_num_steps),
        coupled_system_(f, y0_, msgs, args...),
        coupled_state_(coupled_system_.initial_state()) {
    check_finite(function_name, "initial state", y0_);
    check_finite(function_name, "initial time", t0_);
    check_finite(function_name, "times", ts_);

    math::apply(
        [&](const auto&... args_ref) {
          (check_finite(function_name, "ode parameters and data", args_ref),
           ...);
        },
        args_tuple_);

    check_nonzero_size(function_name, "initial state", y0_);
    check_nonzero_size(function_name, "times", ts_);
    check_sorted(function_name, "times", ts_);
    check_less(function_name, "initial time", t0_, ts_[0]);
    check_positive_finite(function_name, "relative_tolerance",
                          relative_tolerance_);
    check_positive_finite(function_name, "absolute_tolerance",
                          absolute_tolerance_);
    check_positive(function_name, "max_num_steps", max_num_steps_);

    nv_state_ = N_VMake_Serial(coupled_state_.size(), coupled_state_.data(),
                               sundials_context_);
  }

  ~arkode_integrator() { N_VDestroy_Serial(nv_state_); }

  /**
   * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
   * times, { t1, t2, t3, ... } using ARKODE's ERKStep embedded explicit
   * Runge-Kutta solver with the Butcher table given by <code>Table</code>.
   *
   * @return std::vector of Eigen::Matrix of the states of the ODE, one for
   *   each solution time (excluding the initial state)
   */
  std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> operator()() {
    std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> y;
    y.reserve(ts_.size());

    void* arkode_mem = ERKStepCreate(&arkode_integrator::erk_rhs, value_of(t0_),
                                     nv_state_, sundials_context_);
    if (arkode_mem == nullptr) {
      throw std::runtime_error("ERKStepCreate failed to allocate memory");
    }

    try {
      CHECK_ARKODE_CALL(
          ERKStepSetUserData(arkode_mem, reinterpret_cast<void*>(this)));
      CHECK_ARKODE_CALL(ERKStepSetTableNum(arkode_mem, Table));
      CHECK_ARKODE_CALL(ERKStepSStolerances(arkode_mem, relative_tolerance_,
                                            absolute_tolerance_));
      CHECK_ARKODE_CALL(ERKStepSetMaxNumSteps(arkode_mem, max_num_steps_));

      // Cap the internal step size at the smallest gap between requested
      // output times, so that ARK_NORMAL never has to interpolate far past
      // an output time to reach it -- this keeps the Hermite interpolation
      // error at each output time small relative to local step error.
      double min_gap = std::numeric_limits<double>::infinity();
      double t_prev = value_of(t0_);
      for (size_t n = 0; n < ts_.size(); ++n) {
        double t_n = value_of(ts_[n]);
        min_gap = std::min(min_gap, t_n - t_prev);
        t_prev = t_n;
      }
      if (std::isfinite(min_gap) && min_gap > 0) {
        CHECK_ARKODE_CALL(ERKStepSetMaxStep(arkode_mem, min_gap));
      }

      double t_init = value_of(t0_);
      for (size_t n = 0; n < ts_.size(); ++n) {
        double t_final = value_of(ts_[n]);

        if (t_final != t_init) {
          int flag = ERKStepEvolve(arkode_mem, t_final, nv_state_, &t_init,
                                   ARK_NORMAL);
          if (flag == ARK_TOO_MUCH_WORK) {
            throw_domain_error(function_name_, "", t_final,
                               "Failed to integrate to next output time (",
                               ") in less than max_num_steps steps");
          }
          CHECK_ARKODE_CALL(flag);
        }

        y.emplace_back(math::apply(
            [&](const auto&... args_ref) {
              return ode_store_sensitivities(f_, coupled_state_, y0_, t0_,
                                             ts_[n], msgs_, args_ref...);
            },
            args_tuple_));

        t_init = t_final;
      }
    } catch (const std::exception& e) {
      ERKStepFree(&arkode_mem);
      throw;
    }

    ERKStepFree(&arkode_mem);

    return y;
  }
};  // arkode_integrator

}  // namespace math
}  // namespace stan
#endif
