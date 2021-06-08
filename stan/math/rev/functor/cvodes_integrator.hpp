#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/rev/functor/coupled_ode_system.hpp>
#include <stan/math/rev/functor/ode_store_sensitivities.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Integrator interface for CVODES' ODE solvers (Adams & BDF
 * methods).
 *
 * @tparam Lmm ID of ODE solver (1: ADAMS, 2: BDF)
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_param Type of scalars for parameters
 * @tparam T_t0 Type of scalar of initial time point
 * @tparam T_ts Type of time-points where ODE solution is returned
 */
template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator {
  using T_Return = return_type_t<T_y0, T_t0, T_ts, T_Args...>;
  using T_y0_t0 = return_type_t<T_y0, T_t0>;

  const char* function_name_;
  const F& f_;
  const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> y0_;
  const T_t0 t0_;
  const std::vector<T_ts>& ts_;
  std::tuple<const T_Args&...> args_tuple_;
  std::tuple<plain_type_t<decltype(value_of(std::declval<const T_Args&>()))>...>
      value_of_args_tuple_;
  const size_t N_;
  std::ostream* msgs_;
  double relative_tolerance_;
  double absolute_tolerance_;
  long int max_num_steps_;  // NOLINT(runtime/int)

  const size_t num_y0_vars_;
  const size_t num_args_vars_;

  coupled_ode_system<F, T_y0_t0, T_Args...> coupled_ode_;

  std::vector<double> coupled_state_;
  N_Vector nv_state_;
  N_Vector* nv_state_sens_;
  SUNMatrix A_;
  SUNLinearSolver LS_;

  /**
   * Implements the function of type CVRhsFn which is the user-defined
   * ODE RHS passed to CVODES.
   */
  static int cv_rhs(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    cvodes_integrator* integrator = static_cast<cvodes_integrator*>(user_data);
    integrator->rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return 0;
  }

  /**
   * Implements the function of type CVSensRhsFn which is the
   * RHS of the sensitivity ODE system.
   */
  static int cv_rhs_sens(int Ns, realtype t, N_Vector y, N_Vector ydot,
                         N_Vector* yS, N_Vector* ySdot, void* user_data,
                         N_Vector tmp1, N_Vector tmp2) {
    cvodes_integrator* integrator = static_cast<cvodes_integrator*>(user_data);
    integrator->rhs_sens(t, NV_DATA_S(y), yS, ySdot);
    return 0;
  }

  /**
   * Implements the function of type CVDlsJacFn which is the
   * user-defined callback for CVODES to calculate the jacobian of the
   * ode_rhs wrt to the states y. The jacobian is stored in column
   * major format.
   */
  static int cv_jacobian_states(realtype t, N_Vector y, N_Vector fy,
                                SUNMatrix J, void* user_data, N_Vector tmp1,
                                N_Vector tmp2, N_Vector tmp3) {
    cvodes_integrator* integrator = static_cast<cvodes_integrator*>(user_data);
    integrator->jacobian_states(t, NV_DATA_S(y), J);
    return 0;
  }

  /**
   * Calculates the ODE RHS, dy_dt, using the user-supplied functor at
   * the given time t and state y.
   */
  inline void rhs(double t, const double y[], double dy_dt[]) const {
    const Eigen::VectorXd y_vec = Eigen::Map<const Eigen::VectorXd>(y, N_);

    Eigen::VectorXd dy_dt_vec
        = apply([&](auto&&... args) { return f_(t, y_vec, msgs_, args...); },
                value_of_args_tuple_);

    check_size_match("cvodes_integrator", "dy_dt", dy_dt_vec.size(), "states",
                     N_);

    std::copy(dy_dt_vec.data(), dy_dt_vec.data() + dy_dt_vec.size(), dy_dt);
  }

  /**
   * Calculates the jacobian of the ODE RHS wrt to its states y at the
   * given time-point t and state y.
   */
  inline void jacobian_states(double t, const double y[], SUNMatrix J) const {
    Eigen::VectorXd fy;
    Eigen::MatrixXd Jfy;

    auto f_wrapped = [&](const Eigen::Matrix<var, Eigen::Dynamic, 1>& y) {
      return apply([&](auto&&... args) { return f_(t, y, msgs_, args...); },
                   value_of_args_tuple_);
    };

    jacobian(f_wrapped, Eigen::Map<const Eigen::VectorXd>(y, N_), fy, Jfy);

    for (size_t j = 0; j < Jfy.cols(); ++j) {
      for (size_t i = 0; i < Jfy.rows(); ++i) {
        SM_ELEMENT_D(J, i, j) = Jfy(i, j);
      }
    }
  }

  /**
   * Calculates the RHS of the sensitivity ODE system which
   * corresponds to the coupled ode system from which the first N
   * states are omitted, since the first N states are the ODE RHS
   * which CVODES separates from the main ODE RHS.
   */
  inline void rhs_sens(double t, const double y[], N_Vector* yS,
                       N_Vector* ySdot) {
    std::vector<double> z(coupled_state_.size());
    std::vector<double> dz_dt;
    std::copy(y, y + N_, z.data());
    for (std::size_t s = 0; s < num_y0_vars_ + num_args_vars_; s++) {
      std::copy(NV_DATA_S(yS[s]), NV_DATA_S(yS[s]) + N_,
                z.data() + (s + 1) * N_);
    }
    coupled_ode_(z, dz_dt, t);
    for (std::size_t s = 0; s < num_y0_vars_ + num_args_vars_; s++) {
      std::move(dz_dt.data() + (s + 1) * N_, dz_dt.data() + (s + 2) * N_,
                NV_DATA_S(ySdot[s]));
    }
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
   * @param relative_tolerance Relative tolerance passed to CVODES
   * @param absolute_tolerance Absolute tolerance passed to CVODES
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
  cvodes_integrator(const char* function_name, const F& f, const T_y0& y0,
                    const T_t0& t0, const std::vector<T_ts>& ts,
                    double relative_tolerance, double absolute_tolerance,
                    long int max_num_steps,  // NOLINT(runtime/int)
                    std::ostream* msgs, const T_Args&... args)
      : function_name_(function_name),
        f_(f),
        y0_(y0.template cast<T_y0_t0>()),
        t0_(t0),
        ts_(ts),
        args_tuple_(args...),
        value_of_args_tuple_(value_of(args)...),
        N_(y0.size()),
        msgs_(msgs),
        relative_tolerance_(relative_tolerance),
        absolute_tolerance_(absolute_tolerance),
        max_num_steps_(max_num_steps),
        num_y0_vars_(count_vars(y0_)),
        num_args_vars_(count_vars(args...)),
        coupled_ode_(f, y0_, msgs, args...),
        coupled_state_(coupled_ode_.initial_state()) {
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
    check_positive_finite(function_name, "relative_tolerance",
                          relative_tolerance_);
    check_positive_finite(function_name, "absolute_tolerance",
                          absolute_tolerance_);
    check_positive(function_name, "max_num_steps", max_num_steps_);

    nv_state_ = N_VMake_Serial(N_, &coupled_state_[0]);
    nv_state_sens_ = nullptr;
    A_ = SUNDenseMatrix(N_, N_);
    LS_ = SUNDenseLinearSolver(nv_state_, A_);

    if (num_y0_vars_ + num_args_vars_ > 0) {
      nv_state_sens_ = N_VCloneVectorArrayEmpty_Serial(
          num_y0_vars_ + num_args_vars_, nv_state_);
      for (std::size_t i = 0; i < num_y0_vars_ + num_args_vars_; i++) {
        NV_DATA_S(nv_state_sens_[i]) = &coupled_state_[N_] + i * N_;
      }
    }
  }

  ~cvodes_integrator() {
    SUNLinSolFree(LS_);
    SUNMatDestroy(A_);
    N_VDestroy_Serial(nv_state_);
    if (num_y0_vars_ + num_args_vars_ > 0) {
      N_VDestroyVectorArray_Serial(nv_state_sens_,
                                   num_y0_vars_ + num_args_vars_);
    }
  }

  /**
   * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
   * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
   * (BDF) solver in CVODES.
   *
   * @return std::vector of Eigen::Matrix of the states of the ODE, one for each
   *   solution time (excluding the initial state)
   */
  std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> operator()() {
    std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> y;

    void* cvodes_mem = CVodeCreate(Lmm);
    if (cvodes_mem == nullptr) {
      throw std::runtime_error("CVodeCreate failed to allocate memory");
    }

    try {
      check_flag_sundials(CVodeInit(cvodes_mem, &cvodes_integrator::cv_rhs,
                                    value_of(t0_), nv_state_),
                          "CVodeInit");

      // Assign pointer to this as user data
      check_flag_sundials(
          CVodeSetUserData(cvodes_mem, reinterpret_cast<void*>(this)),
          "CVodeSetUserData");

      cvodes_set_options(cvodes_mem, max_num_steps_);

      check_flag_sundials(CVodeSStolerances(cvodes_mem, relative_tolerance_,
                                            absolute_tolerance_),
                          "CVodeSStolerances");

      check_flag_sundials(CVodeSetLinearSolver(cvodes_mem, LS_, A_),
                          "CVodeSetLinearSolver");
      check_flag_sundials(
          CVodeSetJacFn(cvodes_mem, &cvodes_integrator::cv_jacobian_states),
          "CVodeSetJacFn");

      // initialize forward sensitivity system of CVODES as needed
      if (num_y0_vars_ + num_args_vars_ > 0) {
        check_flag_sundials(
            CVodeSensInit(
                cvodes_mem, static_cast<int>(num_y0_vars_ + num_args_vars_),
                CV_STAGGERED, &cvodes_integrator::cv_rhs_sens, nv_state_sens_),
            "CVodeSensInit");

        check_flag_sundials(CVodeSetSensErrCon(cvodes_mem, SUNTRUE),
                            "CVodeSetSensErrCon");

        check_flag_sundials(CVodeSensEEtolerances(cvodes_mem),
                            "CVodeSensEEtolerances");
      }

      double t_init = value_of(t0_);
      for (size_t n = 0; n < ts_.size(); ++n) {
        double t_final = value_of(ts_[n]);

        if (t_final != t_init) {
          int error_code
              = CVode(cvodes_mem, t_final, nv_state_, &t_init, CV_NORMAL);

          if (error_code == CV_TOO_MUCH_WORK) {
            throw_domain_error(function_name_, "", t_final,
                               "Failed to integrate to next output time (",
                               ") in less than max_num_steps steps");
          } else {
            check_flag_sundials(error_code, "CVode");
          }

          if (num_y0_vars_ + num_args_vars_ > 0) {
            check_flag_sundials(
                CVodeGetSens(cvodes_mem, &t_init, nv_state_sens_),
                "CVodeGetSens");
          }
        }

        y.emplace_back(apply(
            [&](auto&&... args) {
              return ode_store_sensitivities(f_, coupled_state_, y0_, t0_,
                                             ts_[n], msgs_, args...);
            },
            args_tuple_));

        t_init = t_final;
      }
    } catch (const std::exception& e) {
      CVodeFree(&cvodes_mem);
      throw;
    }

    CVodeFree(&cvodes_mem);

    return y;
  }
};  // cvodes integrator

}  // namespace math
}  // namespace stan
#endif
