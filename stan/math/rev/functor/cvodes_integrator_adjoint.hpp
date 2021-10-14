#ifndef STAN_MATH_REV_FUNCTOR_CVODES_INTEGRATOR_ADJOINT_HPP
#define STAN_MATH_REV_FUNCTOR_CVODES_INTEGRATOR_ADJOINT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/save_varis.hpp>
#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/rev/functor/ode_store_sensitivities.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Integrator interface for CVODES' adjoint ODE solvers (Adams & BDF
 * methods).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of scalars for initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of time-points where ODE solution is returned
 * @tparam T_Args Types of pass-through parameters
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args>
class cvodes_integrator_adjoint_vari : public vari_base {
  using T_Return = return_type_t<T_y0, T_t0, T_ts, T_Args...>;
  using T_y0_t0 = return_type_t<T_y0, T_t0>;

  static constexpr bool is_var_ts_{is_var<T_ts>::value};
  static constexpr bool is_var_t0_{is_var<T_t0>::value};
  static constexpr bool is_var_y0_{is_var<T_y0>::value};
  static constexpr bool is_var_y0_t0_{is_var<T_y0_t0>::value};
  static constexpr bool is_any_var_args_{
      disjunction<is_var<scalar_type_t<T_Args>>...>::value};
  static constexpr bool is_var_return_{is_var<T_Return>::value};
  static constexpr bool is_var_only_ts_{
      is_var_ts_ && !(is_var_t0_ || is_var_y0_t0_ || is_any_var_args_)};

  arena_t<std::vector<Eigen::VectorXd>> y_;
  arena_t<std::vector<T_ts>> ts_;
  arena_t<Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1>> y0_;
  arena_t<Eigen::VectorXd> absolute_tolerance_backward_;
  arena_t<Eigen::VectorXd> state_backward_;
  size_t num_args_vars_;
  arena_t<Eigen::VectorXd> quad_;
  arena_t<T_t0> t0_;
  arena_t<std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>>> y_return_;

  double relative_tolerance_forward_;
  double relative_tolerance_backward_;
  double relative_tolerance_quadrature_;
  double absolute_tolerance_quadrature_;
  long int max_num_steps_;                  // NOLINT(runtime/int)
  long int num_steps_between_checkpoints_;  // NOLINT(runtime/int)
  size_t N_;
  std::ostream* msgs_;
  vari** args_varis_;
  int interpolation_polynomial_;
  int solver_forward_;
  int solver_backward_;
  int index_backward_;
  bool backward_is_initialized_{false};

  /**
   * Since the CVODES solver manages memory with malloc calls, these resources
   * must be freed using a destructor call (which is not being called for the
   * vari class).
   */
  struct cvodes_solver : public chainable_alloc {
    const char* function_name_;
    std::decay_t<F> f_;
    const size_t N_;
    N_Vector nv_state_backward_;
    N_Vector nv_quad_;
    N_Vector nv_absolute_tolerance_backward_;
    SUNMatrix A_backward_;
    SUNLinearSolver LS_backward_;
    void* cvodes_mem_;
    std::tuple<T_Args...> local_args_tuple_;
    using value_tuple_t = std::tuple<
        promote_scalar_t<partials_type_t<scalar_type_t<T_Args>>, T_Args>...>;
    value_tuple_t value_of_args_tuple_;

    template <typename FF, typename StateBwd, typename Quad, typename AbsTolBwd>
    cvodes_solver(const char* function_name, FF&& f, size_t N,
                  size_t num_args_vars, size_t ts_size, int solver_forward,
                  StateBwd& state_backward, Quad& quad,
                  AbsTolBwd& absolute_tolerance_backward, const T_Args&... args)
        : chainable_alloc(),
          function_name_(function_name),
          f_(std::forward<FF>(f)),
          N_(N),
          nv_state_backward_(N_VMake_Serial(N, state_backward.data())),
          nv_quad_(N_VMake_Serial(num_args_vars, quad.data())),
          nv_absolute_tolerance_backward_(
              N_VMake_Serial(N, absolute_tolerance_backward.data())),
          A_backward_(SUNDenseMatrix(N, N)),
          LS_backward_(
              N == 0 ? nullptr
                     : SUNDenseLinearSolver(nv_state_backward_, A_backward_)),
          cvodes_mem_([](int solver_forward) {
            void* cvodes_mem = CVodeCreate(solver_forward);
            if (cvodes_mem == nullptr) {
              throw std::runtime_error("CVodeCreate failed to allocate memory");
            }
            return cvodes_mem;
          }(solver_forward)),
          local_args_tuple_(deep_copy_vars(args)...),
          value_of_args_tuple_(apply(
              [](auto&&... args) { return value_tuple_t(value_of(args)...); },
              local_args_tuple_)) {}

    virtual ~cvodes_solver() {
      SUNMatDestroy(A_backward_);
      if (N_ > 0) {
        SUNLinSolFree(LS_backward_);
      }
      N_VDestroy_Serial(nv_state_backward_);
      N_VDestroy_Serial(nv_quad_);
      N_VDestroy_Serial(nv_absolute_tolerance_backward_);

      CVodeFree(&cvodes_mem_);
    }
  };
  cvodes_solver* solver_{nullptr};

 public:
  /**
   * Construct cvodes_integrator object. Note: All arguments must be stored as
   * copies if in doubt. The reason is that the references can go out of scope,
   * since the work done from the integrator is in the chain method.
   *
   * @param function_name Calling function name (for printing debugging
   * messages)
   * @param f Right hand side of the ODE
   * @param y0 Initial state
   * @param t0 Initial time
   * @param ts Times at which to solve the ODE at. All values must be sorted and
   *   not less than t0.
   * @param relative_tolerance_forward Relative tolerance for forward problem
   * passed to CVODES
   * @param absolute_tolerance_forward Absolute tolerance per ODE state for
   * forward problem passed to CVODES
   * @param relative_tolerance_backward Relative tolerance for backward problem
   * passed to CVODES
   * @param absolute_tolerance_backward Absolute tolerance per ODE state for
   * backward problem passed to CVODES
   * @param relative_tolerance_quadrature Relative tolerance for quadrature
   * problem passed to CVODES
   * @param absolute_tolerance_quadrature Absolute tolerance for quadrature
   * problem passed to CVODES
   * @param max_num_steps Upper limit on the number of integration steps to
   *   take between each output (error if exceeded)
   * @param num_steps_between_checkpoints Number of integrator steps after which
   * a checkpoint is stored for the backward pass
   * @param interpolation_polynomial type of polynomial used for interpolation
   * @param solver_forward solver used for forward pass
   * @param solver_backward solver used for backward pass
   *   take between each output (error if exceeded)
   * @param[in, out] msgs the print stream for warning messages
   * @param args Extra arguments passed unmodified through to ODE right hand
   * side function
   * @return Solution to ODE at times \p ts
   * @return a vector of states, each state being a vector of the
   * same size as the state variable, corresponding to a time in ts.
   */
  template <typename FF, typename T_y00,
            require_eigen_col_vector_t<T_y00>* = nullptr>
  cvodes_integrator_adjoint_vari(
      const char* function_name, FF&& f, const T_y00& y0, const T_t0& t0,
      const std::vector<T_ts>& ts, double relative_tolerance_forward,
      const Eigen::VectorXd& absolute_tolerance_forward,
      double relative_tolerance_backward,
      const Eigen::VectorXd& absolute_tolerance_backward,
      double relative_tolerance_quadrature,
      double absolute_tolerance_quadrature,
      long int max_num_steps,                  // NOLINT(runtime/int)
      long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
      int interpolation_polynomial, int solver_forward, int solver_backward,
      std::ostream* msgs, const T_Args&... args)
      : vari_base(),
        y_(ts.size()),
        ts_([](const char* function_name, auto&& ts) {
          check_nonzero_size(function_name, "times", ts);
          check_finite(function_name, "times", ts);
          check_sorted(function_name, "times", ts);
          return arena_t<std::vector<T_ts>>(ts.begin(), ts.end());
        }(function_name, ts)),
        y0_([](const char* function_name, auto&& y0) {
          check_nonzero_size(function_name, "initial state", y0);
          check_finite(function_name, "initial state", y0);
          return y0;
        }(function_name, y0)),
        absolute_tolerance_backward_([](const char* function_name, auto N,
                                        auto&& absolute_tolerance_backward) {
          check_positive_finite(function_name, "absolute_tolerance_backward",
                                absolute_tolerance_backward);
          check_size_match(function_name, "absolute_tolerance_backward",
                           absolute_tolerance_backward.size(), "states", N);
          return absolute_tolerance_backward;
        }(function_name, y0_.size(), absolute_tolerance_backward)),
        state_backward_(Eigen::VectorXd::Zero(y0_.size())),
        num_args_vars_(count_vars(args...)),
        quad_(Eigen::VectorXd::Zero(num_args_vars_)),
        t0_([](const char* function_name, auto&& t0, auto&& ts) {
          check_finite(function_name, "initial time", t0);
          check_less(function_name, "initial time", t0, ts[0]);
          return t0;
        }(function_name, t0, ts)),
        y_return_(ts.size()),
        relative_tolerance_forward_(
            [](const char* function_name, auto&& relative_tolerance_forward) {
              check_positive_finite(function_name, "relative_tolerance_forward",
                                    relative_tolerance_forward);
              return relative_tolerance_forward;
            }(function_name, relative_tolerance_forward)),
        relative_tolerance_backward_([](const char* function_name,
                                        auto&& relative_tolerance_backward) {
          check_positive_finite(function_name, "relative_tolerance_backward",
                                relative_tolerance_backward);
          return relative_tolerance_backward;
        }(function_name, relative_tolerance_backward)),
        relative_tolerance_quadrature_(
            [](const char* function_name,
               auto&& relative_tolerance_quadrature) {
              check_positive_finite(function_name,
                                    "relative_tolerance_quadrature",
                                    relative_tolerance_quadrature);
              return relative_tolerance_quadrature;
            }(function_name, relative_tolerance_quadrature)),
        absolute_tolerance_quadrature_(
            [](const char* function_name,
               auto&& absolute_tolerance_quadrature) {
              check_positive_finite(function_name,
                                    "absolute_tolerance_quadrature",
                                    absolute_tolerance_quadrature);
              return absolute_tolerance_quadrature;
            }(function_name, absolute_tolerance_quadrature)),
        max_num_steps_([](const char* function_name, auto&& max_num_steps) {
          check_positive(function_name, "max_num_steps", max_num_steps);
          return max_num_steps;
        }(function_name, max_num_steps)),
        num_steps_between_checkpoints_(
            [](const char* function_name,
               auto&& num_steps_between_checkpoints) {
              check_positive(function_name, "num_steps_between_checkpoints",
                             num_steps_between_checkpoints);
              return num_steps_between_checkpoints;
            }(function_name, num_steps_between_checkpoints)),
        N_(y0_.size()),
        msgs_(msgs),
        args_varis_([](auto num_vars, auto&&... args) {
          vari** vari_mem
              = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
                  num_vars);
          save_varis(vari_mem, args...);
          return vari_mem;
        }(this->num_args_vars_, args...)),
        interpolation_polynomial_([](const char* function_name,
                                     auto&& interpolation_polynomial) {
          // for polynomial: 1=CV_HERMITE / 2=CV_POLYNOMIAL
          if (interpolation_polynomial != 1 && interpolation_polynomial != 2)
            invalid_argument(function_name, "interpolation_polynomial",
                             interpolation_polynomial, "",
                             ", must be 1 for Hermite or 2 for polynomial "
                             "interpolation of ODE solution");
          return interpolation_polynomial;
        }(function_name, interpolation_polynomial)),
        solver_forward_([](const char* function_name, auto&& solver_forward) {
          // 1=Adams=CV_ADAMS, 2=BDF=CV_BDF
          if (solver_forward != 1 && solver_forward != 2)
            invalid_argument(
                function_name, "solver_forward", solver_forward, "",
                ", must be 1 for Adams or 2 for BDF forward solver");
          return solver_forward;
        }(function_name, solver_forward)),
        solver_backward_([](const char* function_name, auto&& solver_backward) {
          if (solver_backward != 1 && solver_backward != 2)
            invalid_argument(
                function_name, "solver_backward", solver_backward, "",
                ", must be 1 for Adams or 2 for BDF backward solver");
          return solver_backward;
        }(function_name, solver_backward)),
        backward_is_initialized_(false),
        solver_(new cvodes_solver(function_name, std::forward<FF>(f), N_,
                                  num_args_vars_, ts_.size(), solver_forward_,
                                  state_backward_, quad_,
                                  absolute_tolerance_backward_, args...)) {
    {
      int param_number = 1;
      stan::math::for_each(
          [&param_number, func_name = function_name](auto&& arg) {
            check_finite(
                func_name,
                (std::string("ode parameters and data argument number ")
                 + std::to_string(param_number) + " at index ")
                    .c_str(),
                arg);
            param_number++;
          },
          solver_->local_args_tuple_);
    }
    Eigen::VectorXd absolute_tolerance_forward_(absolute_tolerance_forward);
    check_positive_finite(function_name, "absolute_tolerance_forward",
                          absolute_tolerance_forward_);
    check_size_match(function_name, "absolute_tolerance_forward",
                     absolute_tolerance_forward_.size(), "states", N_);
    N_Vector nv_absolute_tolerance_forward_;
    SUNLinearSolver LS_forward_;
    SUNMatrix A_forward_;
    N_Vector nv_state_forward_;
    Eigen::VectorXd state_forward_ = value_of(y0_);

    try {
      nv_state_forward_ = N_VMake_Serial(N_, state_forward_.data());
      nv_absolute_tolerance_forward_
          = N_VMake_Serial(N_, absolute_tolerance_forward_.data());
      A_forward_ = SUNDenseMatrix(N_, N_);
      LS_forward_ = N_ == 0
                        ? nullptr
                        : SUNDenseLinearSolver(nv_state_forward_, A_forward_);

      check_flag_sundials(CVodeInit(solver_->cvodes_mem_,
                                    &cvodes_integrator_adjoint_vari::cv_rhs,
                                    value_of(t0_), nv_state_forward_),
                          "CVodeInit");

      // Assign pointer to this as user data
      check_flag_sundials(
          CVodeSetUserData(solver_->cvodes_mem_, reinterpret_cast<void*>(this)),
          "CVodeSetUserData");

      cvodes_set_options(solver_->cvodes_mem_, max_num_steps_);

      check_flag_sundials(
          CVodeSVtolerances(solver_->cvodes_mem_, relative_tolerance_forward_,
                            nv_absolute_tolerance_forward_),
          "CVodeSVtolerances");

      check_flag_sundials(
          CVodeSetLinearSolver(solver_->cvodes_mem_, LS_forward_, A_forward_),
          "CVodeSetLinearSolver");

      check_flag_sundials(
          CVodeSetJacFn(
              solver_->cvodes_mem_,
              &cvodes_integrator_adjoint_vari::cv_jacobian_rhs_states),
          "CVodeSetJacFn");

      // initialize backward sensitivity system of CVODES as needed
      if (is_var_return_ && !is_var_only_ts_) {
        check_flag_sundials(
            CVodeAdjInit(solver_->cvodes_mem_, num_steps_between_checkpoints_,
                         interpolation_polynomial_),
            "CVodeAdjInit");
      }

      /**
       * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set
       * of times, { t1, t2, t3, ... } using the requested forward solver of
       * CVODES.
       */
      const auto ts_dbl = value_of(ts_);

      double t_init = value_of(t0_);
      for (size_t n = 0; n < ts_dbl.size(); ++n) {
        double t_final = ts_dbl[n];
        if (t_final != t_init) {
          if (is_var_return_ && !is_var_only_ts_) {
            int ncheck;

            int error_code
                = CVodeF(solver_->cvodes_mem_, t_final, nv_state_forward_,
                         &t_init, CV_NORMAL, &ncheck);

            if (unlikely(error_code == CV_TOO_MUCH_WORK)) {
              throw_domain_error(solver_->function_name_, "", t_final,
                                 "Failed to integrate to next output time (",
                                 ") in less than max_num_steps steps");
            } else {
              check_flag_sundials(error_code, "CVodeF");
            }
          } else {
            int error_code = CVode(solver_->cvodes_mem_, t_final,
                                   nv_state_forward_, &t_init, CV_NORMAL);

            if (unlikely(error_code == CV_TOO_MUCH_WORK)) {
              throw_domain_error(solver_->function_name_, "", t_final,
                                 "Failed to integrate to next output time (",
                                 ") in less than max_num_steps steps");
            } else {
              check_flag_sundials(error_code, "CVode");
            }
          }
        }
        // Need to make a deep copy of state_forward_
        auto state_forward_tmp = Eigen::VectorXd(state_forward_);
        y_[n] = state_forward_tmp;
        this->y_return_[n] = state_forward_tmp;
        t_init = t_final;
      }
      N_VDestroy_Serial(nv_absolute_tolerance_forward_);
      SUNMatDestroy(A_forward_);
      if (N_ > 0) {
        SUNLinSolFree(LS_forward_);
      }
      N_VDestroy_Serial(nv_state_forward_);
      ChainableStack::instance_->var_stack_.push_back(this);

    } catch (...) {
      N_VDestroy_Serial(nv_absolute_tolerance_forward_);
      SUNMatDestroy(A_forward_);
      if (N_ > 0) {
        SUNLinSolFree(LS_forward_);
      }
      N_VDestroy_Serial(nv_state_forward_);
      std::rethrow_exception(std::current_exception());
    }
  }

 public:
  /**
   * Obtain solution of ODE.
   *
   * @return std::vector of Eigen::Matrix of the states of the ODE, one for each
   *   solution time (excluding the initial state)
   */
  std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> solution() noexcept {
    return std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>>(
        this->y_return_.begin(), this->y_return_.end());
  }

  /**
   * No-op for setting adjoints since this class does not own any adjoints.
   */
  void set_zero_adjoint() final{};

  void chain() final {
    if (!is_var_return_) {
      return;
    }
    using arena_var_vector = arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>>;
    // for sensitivities wrt to ts we do not need to run the backward
    // integration
    if (is_var_ts_) {
      for (int i = 0; i < ts_.size(); ++i) {
        adjoint_of(ts_[i]) += forward_as<arena_var_vector>(this->y_return_[i])
                                  .adj()
                                  .dot(rhs(value_of(ts_[i]), y_[i],
                                           solver_->value_of_args_tuple_));
      }

      if (is_var_only_ts_) {
        return;
      }
    }

    state_backward_.setZero();
    quad_.setZero();

    // At every time step, collect the adjoints from the output
    // variables and re-initialize the solver
    double t_init = value_of(ts_.back());
    for (int i = ts_.size() - 1; i >= 0; --i) {
      // Take in the adjoints from all the output variables at this point
      // in time
      state_backward_ += forward_as<arena_var_vector>(this->y_return_[i]).adj();

      double t_final = value_of((i > 0) ? ts_[i - 1] : t0_);
      if (t_final != t_init) {
        if (unlikely(!backward_is_initialized_)) {
          check_flag_sundials(CVodeCreateB(solver_->cvodes_mem_,
                                           solver_backward_, &index_backward_),
                              "CVodeCreateB");

          check_flag_sundials(
              CVodeSetUserDataB(solver_->cvodes_mem_, index_backward_,
                                reinterpret_cast<void*>(this)),
              "CVodeSetUserDataB");

          // initialize CVODES backward machinery.
          // the states of the backward problem *are* the adjoints
          // of the ode states
          check_flag_sundials(
              CVodeInitB(solver_->cvodes_mem_, index_backward_,
                         &cvodes_integrator_adjoint_vari::cv_rhs_adj, t_init,
                         solver_->nv_state_backward_),
              "CVodeInitB");

          check_flag_sundials(
              CVodeSVtolerancesB(solver_->cvodes_mem_, index_backward_,
                                 relative_tolerance_backward_,
                                 solver_->nv_absolute_tolerance_backward_),
              "CVodeSVtolerancesB");

          check_flag_sundials(
              CVodeSetMaxNumStepsB(solver_->cvodes_mem_, index_backward_,
                                   max_num_steps_),
              "CVodeSetMaxNumStepsB");

          check_flag_sundials(CVodeSetLinearSolverB(
                                  solver_->cvodes_mem_, index_backward_,
                                  solver_->LS_backward_, solver_->A_backward_),
                              "CVodeSetLinearSolverB");

          check_flag_sundials(
              CVodeSetJacFnB(
                  solver_->cvodes_mem_, index_backward_,
                  &cvodes_integrator_adjoint_vari::cv_jacobian_rhs_adj_states),
              "CVodeSetJacFnB");

          // Allocate space for backwards quadrature needed when
          // parameters vary.
          if (is_any_var_args_) {
            check_flag_sundials(
                CVodeQuadInitB(solver_->cvodes_mem_, index_backward_,
                               &cvodes_integrator_adjoint_vari::cv_quad_rhs_adj,
                               solver_->nv_quad_),
                "CVodeQuadInitB");

            check_flag_sundials(
                CVodeQuadSStolerancesB(solver_->cvodes_mem_, index_backward_,
                                       relative_tolerance_quadrature_,
                                       absolute_tolerance_quadrature_),
                "CVodeQuadSStolerancesB");

            check_flag_sundials(CVodeSetQuadErrConB(solver_->cvodes_mem_,
                                                    index_backward_, SUNTRUE),
                                "CVodeSetQuadErrConB");
          }

          backward_is_initialized_ = true;
        } else {
          // just re-initialize the solver
          check_flag_sundials(
              CVodeReInitB(solver_->cvodes_mem_, index_backward_, t_init,
                           solver_->nv_state_backward_),
              "CVodeReInitB");

          if (is_any_var_args_) {
            check_flag_sundials(
                CVodeQuadReInitB(solver_->cvodes_mem_, index_backward_,
                                 solver_->nv_quad_),
                "CVodeQuadReInitB");
          }
        }

        int error_code = CVodeB(solver_->cvodes_mem_, t_final, CV_NORMAL);

        if (unlikely(error_code == CV_TOO_MUCH_WORK)) {
          throw_domain_error(solver_->function_name_, "", t_final,
                             "Failed to integrate backward to output time (",
                             ") in less than max_num_steps steps");
        } else {
          check_flag_sundials(error_code, "CVodeB");
        }

        // obtain adjoint states and update t_init to time point
        // reached of t_final
        check_flag_sundials(CVodeGetB(solver_->cvodes_mem_, index_backward_,
                                      &t_init, solver_->nv_state_backward_),
                            "CVodeGetB");

        if (is_any_var_args_) {
          check_flag_sundials(
              CVodeGetQuadB(solver_->cvodes_mem_, index_backward_, &t_init,
                            solver_->nv_quad_),
              "CVodeGetQuadB");
        }
      }
    }

    if (is_var_t0_) {
      adjoint_of(t0_) += -state_backward_.dot(
          rhs(t_init, value_of(y0_), solver_->value_of_args_tuple_));
    }

    // After integrating all the way back to t0, we finally have the
    // the adjoints we wanted
    // These are the dlog_density / d(initial_conditions[s]) adjoints
    if (is_var_y0_t0_) {
      forward_as<arena_var_vector>(y0_).adj() += state_backward_;
    }

    // These are the dlog_density / d(parameters[s]) adjoints
    if (is_any_var_args_) {
      Eigen::Map<Eigen::Matrix<vari*, -1, 1>>(args_varis_, num_args_vars_).adj()
          += quad_;
    }
  }

 private:
  /**
   * Call the ODE RHS with given tuple.
   */
  template <typename yT, typename... ArgsT>
  constexpr auto rhs(double t, const yT& y,
                     const std::tuple<ArgsT...>& args_tuple) const {
    return apply(
        [&](auto&&... args) { return solver_->f_(t, y, msgs_, args...); },
        args_tuple);
  }

  /**
   * Utility to cast user memory pointer passed in from CVODES to actual typed
   * object pointer.
   */
  constexpr static cvodes_integrator_adjoint_vari* cast_to_self(void* mem) {
    return static_cast<cvodes_integrator_adjoint_vari*>(mem);
  }

  /**
   * Calculates the ODE RHS, dy_dt, using the user-supplied functor at
   * the given time t and state y.
   */
  inline int rhs(double t, const double* y, double*& dy_dt) const {
    const Eigen::VectorXd y_vec = Eigen::Map<const Eigen::VectorXd>(y, N_);
    const Eigen::VectorXd dy_dt_vec
        = rhs(t, y_vec, solver_->value_of_args_tuple_);
    check_size_match(solver_->function_name_, "dy_dt", dy_dt_vec.size(),
                     "states", N_);
    Eigen::Map<Eigen::VectorXd>(dy_dt, N_) = dy_dt_vec;
    return 0;
  }

  /**
   * Implements the function of type CVRhsFn which is the user-defined
   * ODE RHS passed to CVODES.
   */
  constexpr static int cv_rhs(realtype t, N_Vector y, N_Vector ydot,
                              void* user_data) {
    return cast_to_self(user_data)->rhs(t, NV_DATA_S(y), NV_DATA_S(ydot));
  }

  /*
   * Calculate the adjoint sensitivity RHS for varying initial conditions
   * and parameters
   *
   * Equation 2.23 in the cvs_guide.
   *
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yB state of the adjoint ODE system
   * @param[out] yBdot evaluation of adjoint ODE RHS
   */
  inline int rhs_adj(double t, N_Vector y, N_Vector yB, N_Vector yBdot) const {
    const nested_rev_autodiff nested;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y_vars(
        Eigen::Map<const Eigen::VectorXd>(NV_DATA_S(y), N_));
    Eigen::Matrix<var, Eigen::Dynamic, 1> f_y_t_vars
        = rhs(t, y_vars, solver_->value_of_args_tuple_);
    check_size_match(solver_->function_name_, "dy_dt", f_y_t_vars.size(),
                     "states", N_);
    f_y_t_vars.adj() = -Eigen::Map<Eigen::VectorXd>(NV_DATA_S(yB), N_);
    grad();
    Eigen::Map<Eigen::VectorXd>(NV_DATA_S(yBdot), N_) = y_vars.adj();
    return 0;
  }

  /**
   * Implements the function of type CVRhsFnB which is the
   * RHS of the backward ODE system.
   */
  constexpr static int cv_rhs_adj(realtype t, N_Vector y, N_Vector yB,
                                  N_Vector yBdot, void* user_data) {
    return cast_to_self(user_data)->rhs_adj(t, y, yB, yBdot);
  }

  /*
   * Calculate the RHS for the quadrature part of the adjoint ODE
   * problem.
   *
   * This is the integrand of equation 2.22 in the cvs_guide.
   *
   * @param[in] t time
   * @param[in] y state of the base ODE system
   * @param[in] yB state of the adjoint ODE system
   * @param[out] qBdot evaluation of adjoint ODE quadrature RHS
   */
  inline int quad_rhs_adj(double t, N_Vector y, N_Vector yB, N_Vector qBdot) {
    Eigen::Map<const Eigen::VectorXd> y_vec(NV_DATA_S(y), N_);
    const nested_rev_autodiff nested;

    // The vars here do not live on the nested stack so must be zero'd
    // separately

    stan::math::for_each([](auto&& arg) { zero_adjoints(arg); },
                         solver_->local_args_tuple_);
    Eigen::Matrix<var, Eigen::Dynamic, 1> f_y_t_vars
        = rhs(t, y_vec, solver_->local_args_tuple_);
    check_size_match(solver_->function_name_, "dy_dt", f_y_t_vars.size(),
                     "states", N_);
    f_y_t_vars.adj() = -Eigen::Map<Eigen::VectorXd>(NV_DATA_S(yB), N_);
    grad();
    apply(
        [&qBdot](auto&&... args) {
          accumulate_adjoints(NV_DATA_S(qBdot), args...);
        },
        solver_->local_args_tuple_);
    return 0;
  }

  /**
   * Implements the function of type CVQuadRhsFnB which is the
   * RHS of the backward ODE system's quadrature.
   */
  constexpr static int cv_quad_rhs_adj(realtype t, N_Vector y, N_Vector yB,
                                       N_Vector qBdot, void* user_data) {
    return cast_to_self(user_data)->quad_rhs_adj(t, y, yB, qBdot);
  }

  /**
   * Calculates the jacobian of the ODE RHS wrt to its states y at the
   * given time-point t and state y.
   */
  inline int jacobian_rhs_states(double t, N_Vector y, SUNMatrix J) const {
    Eigen::Map<Eigen::MatrixXd> Jfy(SM_DATA_D(J), N_, N_);

    nested_rev_autodiff nested;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y_var(
        Eigen::Map<const Eigen::VectorXd>(NV_DATA_S(y), N_));
    Eigen::Matrix<var, Eigen::Dynamic, 1> fy_var
        = rhs(t, y_var, solver_->value_of_args_tuple_);

    check_size_match(solver_->function_name_, "dy_dt", fy_var.size(), "states",
                     N_);

    grad(fy_var.coeffRef(0).vi_);
    Jfy.col(0) = y_var.adj();
    for (int i = 1; i < fy_var.size(); ++i) {
      nested.set_zero_all_adjoints();
      grad(fy_var.coeffRef(i).vi_);
      Jfy.col(i) = y_var.adj();
    }
    Jfy.transposeInPlace();
    return 0;
  }

  /**
   * Implements the function of type CVDlsJacFn which is the
   * user-defined callback for CVODES to calculate the jacobian of the
   * ode_rhs wrt to the states y. The jacobian is stored in column
   * major format.
   */
  constexpr static int cv_jacobian_rhs_states(realtype t, N_Vector y,
                                              N_Vector fy, SUNMatrix J,
                                              void* user_data, N_Vector tmp1,
                                              N_Vector tmp2, N_Vector tmp3) {
    return cast_to_self(user_data)->jacobian_rhs_states(t, y, J);
  }

  /*
   * Calculate the Jacobian of the RHS of the adjoint ODE (see rhs_adj
   * below for citation for how this is done)
   *
   * @param[in] t Time
   * @param[in] y State of system
   * @param[out] J CVode structure where output is to be stored
   */
  inline int jacobian_rhs_adj_states(double t, N_Vector y, SUNMatrix J) const {
    // J_adj_y = -1 * transpose(J_y)
    int error_code = jacobian_rhs_states(t, y, J);

    Eigen::Map<Eigen::MatrixXd> J_adj_y(SM_DATA_D(J), N_, N_);
    J_adj_y.transposeInPlace();
    J_adj_y.array() *= -1.0;
    return error_code;
  }

  /**
   * Implements the CVLsJacFnB function for evaluating the jacobian of
   * the adjoint problem wrt to the backward states.
   */
  constexpr static int cv_jacobian_rhs_adj_states(realtype t, N_Vector y,
                                                  N_Vector yB, N_Vector fyB,
                                                  SUNMatrix J, void* user_data,
                                                  N_Vector tmp1, N_Vector tmp2,
                                                  N_Vector tmp3) {
    return cast_to_self(user_data)->jacobian_rhs_adj_states(t, y, J);
  }
};  // cvodes integrator adjoint vari

}  // namespace math
}  // namespace stan
#endif
