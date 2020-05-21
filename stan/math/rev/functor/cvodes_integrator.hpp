#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/save_varis.hpp>
#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/rev/functor/coupled_ode_system.hpp>
#include <stan/math/rev/functor/ode_store_sensitivities.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cvodes/cvodes.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename T_Return>
inline std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> build_varis(
    vari* ode_vari, vari**& non_chaining_varis,
    const std::vector<Eigen::VectorXd>& y);

template <>
inline std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> build_varis<var>(
    vari* ode_vari, vari**& non_chaining_varis,
    const std::vector<Eigen::VectorXd>& y) {
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y_return(y.size());

  if (y.size() == 0) {
    return y_return;
  }

  int N = y[0].size();

  non_chaining_varis
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(y.size() * N);

  for (size_t i = 0; i < y.size(); ++i) {
    for (size_t j = 0; j < N; j++) {
      non_chaining_varis[i * N + j] = new vari(y[i][j], false);
    }
  }

  for (size_t i = 0; i < y.size(); ++i) {
    y_return[i].resize(N);
    for (size_t j = 0; j < N; j++)
      y_return[i][j] = var(non_chaining_varis[i * N + j]);
  }

  return y_return;
}

/*
 * If theta and y are both doubles, just pass the values through (there's
 * no autodiff to handle here).
 */
template <>
inline std::vector<Eigen::VectorXd> build_varis<double>(
    vari* ode_vari, vari**& non_chaining_varis,
    const std::vector<Eigen::VectorXd>& y) {
  return y;
}

template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_vari;

template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_memory : public chainable_alloc {
  const size_t N_;
  const F f_;
  const Eigen::Matrix<T_y0, Eigen::Dynamic, 1> y0_;
  const T_t0 t0_;
  const std::vector<T_ts> ts_;
  std::tuple<T_Args...> args_tuple_;
  std::tuple<decltype(value_of(T_Args()))...> value_of_args_tuple_;
  std::vector<Eigen::VectorXd> y_;
  void* cvodes_mem_;
  Eigen::VectorXd state;
  N_Vector nv_state_;
  SUNMatrix A_;
  SUNLinearSolver LS_;

  cvodes_integrator_memory(const F& f,
			   const Eigen::Matrix<T_y0, Eigen::Dynamic, 1>& y0,
			   const T_t0& t0, const std::vector<T_ts>& ts,
			   const T_Args&... args) :
    N_(y0.size()),
    f_(f),
    y0_(y0),
    t0_(t0),
    ts_(ts),
    args_tuple_(std::make_tuple(args...)),
    value_of_args_tuple_(std::make_tuple(value_of(args)...)),
    y_(ts_.size()),
    cvodes_mem_(nullptr),
    state(value_of(y0)) {
    if(N_ > 0) {
      nv_state_ = N_VMake_Serial(N_, state.data());
      A_ = SUNDenseMatrix(N_, N_);
      LS_ = SUNDenseLinearSolver(nv_state_, A_);

      cvodes_mem_ = CVodeCreate(Lmm);
      if (cvodes_mem_ == nullptr) {
	throw std::runtime_error("CVodeCreate failed to allocate memory");
      }
    }
  }

  ~cvodes_integrator_memory() {
    if(N_ > 0) {
      SUNLinSolFree(LS_);
      SUNMatDestroy(A_);

      N_VDestroy_Serial(nv_state_);

      if(cvodes_mem_) {
	CVodeFree(&cvodes_mem_);
      }
    }
  }

  friend class cvodes_integrator_vari<Lmm, F, T_y0, T_t0, T_ts, T_Args...>;
};
  
/**
 * Integrator interface for CVODES' ODE solvers (Adams & BDF
 * methods).
 *
 * @tparam Lmm ID of ODE solver (1: ADAMS, 2: BDF)
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of scalars for initial state
 * @tparam T_param Type of scalars for parameters
 * @tparam T_t0 Type of scalar of initial time point
 * @tparam T_ts Type of time-points where ODE solution is returned
 */
template <int Lmm, typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_vari : public vari {
  using T_Return = return_type_t<T_y0, T_t0, T_ts, T_Args...>;

  const size_t N_;
  bool returned_;
  std::ostream* msgs_;
  double relative_tolerance_;
  double absolute_tolerance_;
  long int max_num_steps_;

  const size_t t0_vars_;
  const size_t ts_vars_;
  const size_t y0_vars_;
  const size_t args_vars_;

  vari** non_chaining_varis_;

  vari** t0_varis_;
  vari** ts_varis_;
  vari** y0_varis_;
  vari** args_varis_;

  cvodes_integrator_memory<Lmm, F, T_y0, T_t0, T_ts, T_Args...>* memory;

  /**
   * Implements the function of type CVRhsFn which is the user-defined
   * ODE RHS passed to CVODES.
   */
  static int cv_rhs(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    const cvodes_integrator_vari* integrator
        = static_cast<const cvodes_integrator_vari*>(user_data);
    integrator->rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
    return 0;
  }

  /**
   * Implements the function of type CVRhsFnB which is the
   * RHS of the backward ODE system.
   */
  static int cv_rhs_adj_sens(realtype t, N_Vector y, N_Vector yB,
                             N_Vector yBdot, void* user_data) {
    const cvodes_integrator_vari* integrator
        = static_cast<const cvodes_integrator_vari*>(user_data);
    integrator->rhs_adj_sens(t, y, yB, yBdot);
    return 0;
  }

  /**
   * Implements the function of type CVQuadRhsFnB which is the
   * RHS of the backward ODE system's quadrature.
   */
  static int cv_quad_rhs_adj(realtype t, N_Vector y, N_Vector yB,
                             N_Vector qBdot, void* user_data) {
    const cvodes_integrator_vari* integrator
        = static_cast<const cvodes_integrator_vari*>(user_data);
    integrator->quad_rhs_adj(t, y, yB, qBdot);
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
    const cvodes_integrator_vari* integrator
        = static_cast<const cvodes_integrator_vari*>(user_data);
    integrator->jacobian_states(t, y, J);
    return 0;
  }

  /**
   * Implements the CVLsJacFnB function for evaluating the jacobian of
   * the adjoint problem.
   */
  static int cv_jacobian_adj(realtype t, N_Vector y, N_Vector yB, N_Vector fyB,
                             SUNMatrix J, void* user_data, N_Vector tmp1,
                             N_Vector tmp2, N_Vector tmp3) {
    const cvodes_integrator_vari* integrator
        = static_cast<const cvodes_integrator_vari*>(user_data);
    integrator->jacobian_adj(t, y, J);
    return 0;
  }

  /**
   * Calculates the ODE RHS, dy_dt, using the user-supplied functor at
   * the given time t and state y.
   */
  inline void rhs(double t, const double y[], double dy_dt[]) const {
    const Eigen::VectorXd y_vec = Eigen::Map<const Eigen::VectorXd>(y, N_);

    Eigen::VectorXd dy_dt_vec
        = apply([&](auto&&... args) { return memory->f_(t, y_vec, msgs_, args...); },
                memory->value_of_args_tuple_);

    check_size_match("cvodes_integrator::rhs", "dy_dt", dy_dt_vec.size(),
                     "states", N_);

    std::copy(dy_dt_vec.data(), dy_dt_vec.data() + dy_dt_vec.size(), dy_dt);
  }

  /*
   * Calculate the adjoint sensitivity RHS for varying initial conditions
   * and parameters
   *
   * @param[in] initial var vector
   * @param[in] param var vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yB state of the adjoint ODE system
   * @param[out] yBdot evaluation of adjoint ODE RHS
   */
  void rhs_adj_sens(double t, N_Vector y, N_Vector yB, N_Vector yBdot) const {
    Eigen::Map<Eigen::VectorXd> y_vec(NV_DATA_S(y), N_);
    Eigen::Map<Eigen::VectorXd> lambda(NV_DATA_S(yB), N_);
    Eigen::Map<Eigen::VectorXd> lambda_dot(NV_DATA_S(yBdot), N_);
    lambda_dot = Eigen::VectorXd::Zero(N_);

    const nested_rev_autodiff nested;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y_vars(y_vec.size());
    for (size_t i = 0; i < y_vars.size(); ++i)
      y_vars(i) = new vari(y_vec(i));

    Eigen::Matrix<var, Eigen::Dynamic, 1> f_y_t_vars
        = apply([&](auto&&... args) { return memory->f_(t, y_vars, msgs_, args...); },
                memory->value_of_args_tuple_);

    check_size_match("coupled_ode_system1", "dy_dt", f_y_t_vars.size(),
                     "states", N_);

    for (size_t i = 0; i < f_y_t_vars.size(); ++i) {
      f_y_t_vars(i).vi_->adj_ = -lambda(i);
    }

    grad();

    for (size_t i = 0; i < y_vars.size(); ++i) {
      lambda_dot(i) = y_vars(i).adj();
    }
  }

  /*
   * Calculate the RHS for the quadrature part of the adjoint ODE problem.
   *
   * @param[in] initial var vector
   * @param[in] param var vector
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yB state of the adjoint ODE system
   * @param[out] qBdot evaluation of adjoint ODE quadrature RHS
   */
  void quad_rhs_adj(double t, N_Vector y, N_Vector yB, N_Vector qBdot) const {
    Eigen::VectorXd y_vec = Eigen::Map<Eigen::VectorXd>(NV_DATA_S(y), N_);
    Eigen::Map<Eigen::VectorXd> lambda(NV_DATA_S(yB), N_);
    Eigen::Map<Eigen::VectorXd> mu_dot(NV_DATA_S(qBdot), args_vars_);
    mu_dot = Eigen::VectorXd::Zero(args_vars_);

    nested_rev_autodiff nested;

    auto local_args_tuple = apply(
				  [&](auto&&... args) {
          return std::tuple<decltype(deep_copy_vars(args))...>(
              deep_copy_vars(args)...);
        },
        memory->args_tuple_);

    Eigen::Matrix<var, Eigen::Dynamic, 1> f_y_t_vars
        = apply([&](auto&&... args) { return memory->f_(t, y_vec, msgs_, args...); },
                local_args_tuple);

    check_size_match("coupled_ode_system2", "dy_dt", f_y_t_vars.size(),
                     "states", N_);

    for (size_t i = 0; i < f_y_t_vars.size(); ++i) {
      f_y_t_vars(i).vi_->adj_ = -lambda(i);
    }

    grad();

    apply([&](auto&&... args) { accumulate_adjoints(mu_dot.data(), args...); },
          local_args_tuple);
  }

  /**
   * Calculates the jacobian of the ODE RHS wrt to its states y at the
   * given time-point t and state y.
   */
  inline void jacobian_states(double t, N_Vector y, SUNMatrix J) const {
    Eigen::VectorXd fy;
    Eigen::MatrixXd Jfy;

    auto f_wrapped = [&](const Eigen::Matrix<var, Eigen::Dynamic, 1>& y) {
      return apply([&](auto&&... args) { return memory->f_(t, y, msgs_, args...); },
                   memory->value_of_args_tuple_);
    };

    jacobian(f_wrapped, Eigen::Map<const Eigen::VectorXd>(NV_DATA_S(y), N_), fy,
             Jfy);

    for (size_t j = 0; j < Jfy.cols(); ++j) {
      for (size_t i = 0; i < Jfy.rows(); ++i) {
        SM_ELEMENT_D(J, i, j) = Jfy(i, j);
      }
    }
  }

  /*
   * Calculate the Jacobian of the RHS of the adjoint ODE (see rhs_adj_sens
   * below for citation for how this is done)
   *
   * @param[in] y State of system
   * @param[in] t Time
   * @param[out] J CVode structure where output is to be stored
   */
  inline void jacobian_adj(double t, N_Vector y, SUNMatrix J) const {
    Eigen::VectorXd fy;
    Eigen::MatrixXd Jfy;

    auto f_wrapped = [&](const Eigen::Matrix<var, Eigen::Dynamic, 1>& y) {
      return apply([&](auto&&... args) { return memory->f_(t, y, msgs_, args...); },
                   memory->value_of_args_tuple_);
    };

    jacobian(f_wrapped, Eigen::Map<const Eigen::VectorXd>(NV_DATA_S(y), N_), fy,
             Jfy);

    for (size_t j = 0; j < Jfy.cols(); ++j) {
      for (size_t i = 0; i < Jfy.rows(); ++i) {
        SM_ELEMENT_D(J, j, i) = -Jfy(i, j);
      }
    }
  }

 public:
  /**
   * Construct cvodes_integrator object
   *
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
   * same size as the state variable, corresponding to a time in ts.
   */
  cvodes_integrator_vari(const F& f,
                         const Eigen::Matrix<T_y0, Eigen::Dynamic, 1>& y0,
                         const T_t0& t0, const std::vector<T_ts>& ts,
                         double relative_tolerance, double absolute_tolerance,
                         long int max_num_steps, std::ostream* msgs,
                         const T_Args&... args)
      : vari(NOT_A_NUMBER),
        N_(y0.size()),
	returned_(false),
	memory(NULL),
        msgs_(msgs),
        relative_tolerance_(relative_tolerance),
        absolute_tolerance_(absolute_tolerance),
        max_num_steps_(max_num_steps),
        t0_vars_(count_vars(t0)),
        ts_vars_(count_vars(ts)),
        y0_vars_(count_vars(y0)),
        args_vars_(count_vars(args...)),
        t0_varis_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(t0_vars_)),
        ts_varis_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(ts_vars_)),
        y0_varis_(
            ChainableStack::instance_->memalloc_.alloc_array<vari*>(y0_vars_)),
        args_varis_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            args_vars_)) {
    const char* fun = "cvodes_integrator::integrate";

    memory = new cvodes_integrator_memory
      <Lmm, F, T_y0, T_t0, T_ts, T_Args...>(f, y0, t0, ts, args...);
  
    save_varis(t0_varis_, t0);
    save_varis(ts_varis_, ts);
    save_varis(y0_varis_, y0);
    save_varis(args_varis_, args...);

    check_finite(fun, "initial state", y0);
    check_finite(fun, "initial time", t0);
    check_finite(fun, "times", ts);

    // Code from: https://stackoverflow.com/a/17340003 . Should probably do
    // something better
    apply(
        [&](auto&&... args) {
          std::vector<int> unused_temp{
              0, (check_finite(fun, "ode parameters and data", args), 0)...};
        },
        memory->args_tuple_);

    check_nonzero_size(fun, "times", ts);
    check_nonzero_size(fun, "initial state", y0);
    check_ordered(fun, "times", ts);
    check_less(fun, "initial time", t0, ts[0]);
    check_positive_finite(fun, "relative_tolerance", relative_tolerance_);
    check_positive_finite(fun, "absolute_tolerance", absolute_tolerance_);
    check_positive(fun, "max_num_steps", max_num_steps_);
  }

  ~cvodes_integrator_vari() {
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
    const double t0_dbl = value_of(memory->t0_);
    const std::vector<double> ts_dbl = value_of(memory->ts_);

    check_flag_sundials(CVodeInit(memory->cvodes_mem_, &cvodes_integrator_vari::cv_rhs, t0_dbl,
				  memory->nv_state_),
			"CVodeInit");

    // Assign pointer to this as user data
    check_flag_sundials(
			CVodeSetUserData(memory->cvodes_mem_, reinterpret_cast<void*>(this)),
			"CVodeSetUserData");

    cvodes_set_options(memory->cvodes_mem_, relative_tolerance_, absolute_tolerance_,
		       max_num_steps_);

    // for the stiff solvers we need to reserve additional memory
    // and provide a Jacobian function call. new API since 3.0.0:
    // create matrix object and linear solver object; resource
    // (de-)allocation is handled in the cvodes_ode_data
    check_flag_sundials(CVodeSetLinearSolver(memory->cvodes_mem_, memory->LS_, memory->A_),
			"CVodeSetLinearSolver");

    check_flag_sundials(
			CVodeSetJacFn(memory->cvodes_mem_,
				      &cvodes_integrator_vari::cv_jacobian_states),
			"CVodeSetJacFn");

    // initialize forward sensitivity system of CVODES as needed
    if (t0_vars_ + ts_vars_ + y0_vars_ + args_vars_ > 0) {
      check_flag_sundials(CVodeAdjInit(memory->cvodes_mem_, 25, CV_HERMITE),
			  "CVodeAdjInit");
    }

    double t_init = t0_dbl;
    for (size_t n = 0; n < ts_dbl.size(); ++n) {
      double t_final = ts_dbl[n];

      if (t_final != t_init) {
	if (t0_vars_ + ts_vars_ + y0_vars_ + args_vars_ > 0) {
	  int ncheck;
	  check_flag_sundials(CVodeF(memory->cvodes_mem_, t_final, memory->nv_state_, &t_init,
				     CV_NORMAL, &ncheck),
			      "CVodeF");
	} else {
	  check_flag_sundials(
			      CVode(memory->cvodes_mem_, t_final, memory->nv_state_, &t_init, CV_NORMAL),
			      "CVode");
	}
      }

      memory->y_[n] = memory->state;

      t_init = t_final;
    }

    returned_ = true;
    return build_varis<T_Return>(this, non_chaining_varis_, memory->y_);
  }

  virtual void chain() {
    // std::cout << "chain" << std::endl; <-- Good way to verify it's only
    //  being called once
    if(memory == NULL)
      return;

    if(memory->cvodes_mem_ == NULL)
      return;
    
    if(returned_ == false)
      return;

    if (t0_vars_ + ts_vars_ + y0_vars_ + args_vars_ == 0) {
      return;
    }

    Eigen::VectorXd state_sens(N_);
    Eigen::VectorXd quad(args_vars_);
    N_Vector nv_state_sens
        = N_VMake_Serial(state_sens.size(), state_sens.data());
    N_Vector nv_quad = N_VMake_Serial(quad.size(), quad.data());
    N_VConst(0.0, nv_state_sens);
    N_VConst(0.0, nv_quad);

    SUNMatrix AB_ = SUNDenseMatrix(N_, N_);
    SUNLinearSolver LSB_ = SUNDenseLinearSolver(nv_state_sens, AB_);

    try {
      int indexB;

      // This is all boilerplate CVODES setting up the adjoint ODE to solve
      check_flag_sundials(CVodeCreateB(memory->cvodes_mem_, Lmm, &indexB),
                          "CVodeCreateB");

      check_flag_sundials(
          CVodeSetUserDataB(memory->cvodes_mem_, indexB, reinterpret_cast<void*>(this)),
          "CVodeSetUserDataB");

      // The ode_rhs_adj_sense functions passed in here cause problems with
      // the autodiff stack (they can cause reallocations of the internal
      // vectors and cause segfaults)
      check_flag_sundials(CVodeInitB(memory->cvodes_mem_, indexB,
                                     &cvodes_integrator_vari::cv_rhs_adj_sens,
                                     value_of(memory->ts_.back()), nv_state_sens),
                          "CVodeInitB");

      check_flag_sundials(
          CVodeSStolerancesB(memory->cvodes_mem_, indexB, relative_tolerance_,
                             absolute_tolerance_),
          "CVodeSStolerancesB");

      check_flag_sundials(CVodeSetLinearSolverB(memory->cvodes_mem_, indexB, LSB_, AB_),
                          "CVodeSetLinearSolverB");

      // The same autodiff issue that applies to ode_rhs_adj_sense applies
      // here
      check_flag_sundials(
          CVodeSetJacFnB(memory->cvodes_mem_, indexB,
                         &cvodes_integrator_vari::cv_jacobian_adj),
          "CVodeSetJacFnB");

      // Allocate space for backwards quadrature
      if (args_vars_ > 0) {
        check_flag_sundials(
            CVodeQuadInitB(memory->cvodes_mem_, indexB,
                           &cvodes_integrator_vari::cv_quad_rhs_adj, nv_quad),
            "CVodeQuadInitB");

        check_flag_sundials(
            CVodeQuadSStolerancesB(memory->cvodes_mem_, indexB, relative_tolerance_,
                                   absolute_tolerance_),
            "CVodeQuadSStolerancesB");

        check_flag_sundials(CVodeSetQuadErrConB(memory->cvodes_mem_, indexB, SUNTRUE),
                            "CVodeSetQuadErrConB");
      }

      // At every time step, collect the adjoints from the output
      // variables and re-initialize the solver
      double t_init = value_of(memory->ts_.back());
      for (int i = memory->ts_.size() - 1; i >= 0; --i) {
        // Take in the adjoints from all the output variables at this point
        // in time
	Eigen::VectorXd step_sens = Eigen::VectorXd::Zero(N_);
        for (int j = 0; j < N_; j++) {
	  //std::cout << "i: " << i << ", j: " << j << std::endl;
          state_sens(j) += non_chaining_varis_[i * N_ + j]->adj_;
          step_sens(j) += non_chaining_varis_[i * N_ + j]->adj_;
        }

        if (ts_vars_ > 0 && i >= 0) {
          ts_varis_[i]->adj_ += apply(
              [&](auto&&... args) {
		double adj = step_sens.dot(memory->f_(t_init, memory->y_[i], msgs_, args...));
		//std::cout << "adj: " << adj << ", i: " << i << std::endl;
                return adj;
              },
              memory->value_of_args_tuple_);
        }

        double t_final = value_of((i > 0) ? memory->ts_[i - 1] : memory->t0_);
        if (t_final != t_init) {
          check_flag_sundials(
              CVodeReInitB(memory->cvodes_mem_, indexB, t_init, nv_state_sens),
              "CVodeReInitB");

	  if(args_vars_ > 0) {
	    check_flag_sundials(CVodeQuadReInitB(memory->cvodes_mem_, indexB, nv_quad),
	    			"CVodeQuadReInitB");
	  }

          check_flag_sundials(CVodeB(memory->cvodes_mem_, t_final, CV_NORMAL),
                              "CVodeB");

	  check_flag_sundials(
	  		      CVodeGetB(memory->cvodes_mem_, indexB, &t_init, nv_state_sens),
	  		      "CVodeGetB");

	  if(args_vars_ > 0) {
	    check_flag_sundials(CVodeGetQuadB(memory->cvodes_mem_, indexB, &t_init, nv_quad),
	    		"CVodeGetQuadB");
	  }
        }
      }

      if (t0_vars_ > 0) {
	Eigen::VectorXd y0d = value_of(memory->y0_);
        t0_varis_[0]->adj_ += apply(
            [&](auto&&... args) {
              return -state_sens.dot(memory->f_(t_init, y0d, msgs_, args...));
            },
            memory->value_of_args_tuple_);
      }

      if (args_vars_ > 0) {
        check_flag_sundials(
            CVodeGetQuadB(memory->cvodes_mem_, indexB, &t_init, nv_quad),
            "CVodeGetQuadB");
      }

      // After integrating all the way back to t0, we finally have the
      // the adjoints we wanted
      // These are the dlog_density / d(initial_conditions[s]) adjoints
      for (size_t s = 0; s < y0_vars_; s++) {
        y0_varis_[s]->adj_ += state_sens(s);
      }

      // These are the dlog_density / d(parameters[s]) adjoints
      for (size_t s = 0; s < args_vars_; s++) {
        args_varis_[s]->adj_ += quad(s);
      }
    } catch (const std::exception& e) {
      SUNLinSolFree(LSB_);
      SUNMatDestroy(AB_);
      N_VDestroy_Serial(nv_state_sens);
      N_VDestroy_Serial(nv_quad);
      throw;
    }

    SUNLinSolFree(LSB_);
    SUNMatDestroy(AB_);
    N_VDestroy_Serial(nv_state_sens);
    N_VDestroy_Serial(nv_quad);
  }
};  // cvodes integrator

}  // namespace math
}  // namespace stan
#endif
