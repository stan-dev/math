#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_ADJOINT_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_ADJOINT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/functor/cvodes_integrator_adjoint.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

  /**
   * Integrator interface for CVODES' ODE solvers (Adams & BDF
   * methods).
   *
   * @tparam F Type of ODE right hand side
   * @tparam T_y0 Type of scalars for initial state
   * @tparam T_t0 Type of initial time
   * @tparam T_ts Type of time-points where ODE solution is returned
   * @tparam T_Args Types of pass-through parameters
   */
template <typename Enable, typename ReturnType, typename F, typename T_y0, typename T_t0, typename T_ts,
  typename... T_Args>
   class cvodes_integrator_adjoint_impl;

   template <typename Enable, typename ReturnType, typename F, typename T_y0_t0, typename T_t0, typename T_ts,
               typename... T_Args>
   struct cvodes_integrator_adjoint_data;


   template <typename ReturnType, typename T_y0, typename T_t0, typename T_ts,
               typename... T_Args>
   struct cvodes_integrator_adjoint_data<require_arithmetic_t<ReturnType>, ReturnType,
     T_y0, T_t0, T_ts, T_Args...> {

     using T_y0_t0 = return_type_t<T_y0, T_t0>;
     std::tuple<arena_t<T_Args>...> local_args_tuple_;
     std::tuple<arena_t<
         promote_scalar_t<partials_type_t<scalar_type_t<T_Args>>, T_Args>>...>
         value_of_args_tuple_;
     /**
      * There's absolutely a way to reuse std::conditional_t<> to figure out
      * based on the templates whether these should be arena memory or not.
      * But we could also just have one for prim and one for rev for now.
      */
     std::vector<Eigen::VectorXd>> y_;
     std::vector<T_ts> ts_;
     Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> y0_;
     Eigen::VectorXd absolute_tolerance_forward_;
     Eigen::VectorXd absolute_tolerance_backward_;
     Eigen::VectorXd state_forward_;
     Eigen::VectorXd state_backward_;
     size_t num_args_vars_;
     Eigen::VectorXd quad_;
     T_t0 t0_;

     double relative_tolerance_forward_;
     double relative_tolerance_backward_;
     double relative_tolerance_quadrature_;
     double absolute_tolerance_quadrature_;
     long int max_num_steps_;                  // NOLINT(runtime/int)
     long int num_steps_between_checkpoints_;  // NOLINT(runtime/int)
     size_t N_;
     std::ostream* msgs_;
     int interpolation_polynomial_;
     int solver_forward_;
     int solver_backward_;
     int index_backward_;
     bool backward_is_initialized_{false};
     template <typename FF, typename T_y0, require_eigen_col_vector_t<T_y0>* = nullptr>
     cvodes_integrator_adjoint_data(
         const char* function_name, FF&& f, const T_y0& y0, const T_t0& t0,
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
         : local_args_tuple_(to_arena(deep_copy_vars(args))...),
           value_of_args_tuple_(to_arena(value_of(args))...),
           y_(ts.size()),
           ts_(ts.begin(), ts.end()),
           y0_(y0),
           absolute_tolerance_forward_(absolute_tolerance_forward),
           absolute_tolerance_backward_(absolute_tolerance_backward),
           state_forward_(value_of(y0)),
           state_backward_(y0.size()),
           num_args_vars_(count_vars(args...)),
           quad_(num_args_vars_),
           t0_(t0),
           relative_tolerance_forward_(relative_tolerance_forward),
           relative_tolerance_backward_(relative_tolerance_backward),
           relative_tolerance_quadrature_(relative_tolerance_quadrature),
           absolute_tolerance_quadrature_(absolute_tolerance_quadrature),
           max_num_steps_(max_num_steps),
           num_steps_between_checkpoints_(num_steps_between_checkpoints),
           N_(y0.size()),
           msgs_(msgs),
           interpolation_polynomial_(interpolation_polynomial),
           solver_forward_(solver_forward),
           solver_backward_(solver_backward),
           backward_is_initialized_(false) {
             check_finite(function_name, "initial state", y0);
             check_finite(function_name, "initial time", t0);
             check_finite(function_name, "times", ts);

             stan::math::for_each(
                 [func_name = function_name](auto&& arg) {
                   check_finite(func_name, "ode parameters and data", arg);
                 },
                 local_args_tuple_);

             check_nonzero_size(function_name, "times", ts);
             check_nonzero_size(function_name, "initial state", y0);
             check_sorted(function_name, "times", ts);
             check_less(function_name, "initial time", t0, ts[0]);
             check_positive_finite(function_name, "relative_tolerance_forward",
                                   relative_tolerance_forward_);
             check_positive_finite(function_name, "absolute_tolerance_forward",
                                   absolute_tolerance_forward_);
             check_size_match(function_name, "absolute_tolerance_forward",
                              absolute_tolerance_forward_.size(), "states", N_);
             check_positive_finite(function_name, "relative_tolerance_backward",
                                   relative_tolerance_backward_);
             check_positive_finite(function_name, "absolute_tolerance_backward",
                                   absolute_tolerance_backward_);
             check_size_match(function_name, "absolute_tolerance_backward",
                              absolute_tolerance_backward_.size(), "states", N_);
             check_positive_finite(function_name, "relative_tolerance_quadrature",
                                   relative_tolerance_quadrature_);
             check_positive_finite(function_name, "absolute_tolerance_quadrature",
                                   absolute_tolerance_quadrature_);
             check_positive(function_name, "max_num_steps", max_num_steps_);
             check_positive(function_name, "num_steps_between_checkpoints",
                            num_steps_between_checkpoints_);
             // for polynomial: 1=CV_HERMITE / 2=CV_POLYNOMIAL
             check_bounded(function_name, "interpolation_polynomial",
                           interpolation_polynomial_, 1, 2);
             // 1=Adams=CV_ADAMS, 2=BDF=CV_BDF
             check_bounded(function_name, "solver_forward", solver_forward_, 1, 2);
             check_bounded(function_name, "solver_backward", solver_backward_, 1, 2);

           }
   };

   /**
    * Overload for all arithmetic types
    * Since the CVODES solver manages memory with malloc calls, these resources
    * must be freed using a destructor call (which is not being called for the
    * vari class).
    */
    template <typename ReturnType>
    struct cvodes_solver {
     std::vector<Eigen::Matrix<ReturnType, Eigen::Dynamic, 1>> y_return_;
     N_Vector nv_state_forward_;
     N_Vector nv_state_backward_;
     N_Vector nv_quad_;
     N_Vector nv_absolute_tolerance_forward_;
     N_Vector nv_absolute_tolerance_backward_;
     SUNMatrix A_forward_;
     SUNLinearSolver LS_forward_;
     SUNMatrix A_backward_;
     SUNLinearSolver LS_backward_;
     const size_t N_;
     const std::string function_name_str_;
     void* cvodes_mem_;
     const std::decay_t<F> f_;

     template <typename FF, typename StateFwd, typename StateBwd, typename Quad,
               typename AbsTolFwd, typename AbsTolBwd>
     cvodes_solver(const char* function_name, FF&& f, size_t N,
                   size_t num_args_vars, size_t ts_size, int solver_forward,
                   StateFwd& state_forward, StateBwd& state_backward, Quad& quad,
                   AbsTolFwd& absolute_tolerance_forward,
                   AbsTolBwd& absolute_tolerance_backward)
         : chainable_alloc(),
           y_return_(ts_size),
           nv_state_forward_(N_VMake_Serial(N, state_forward.data())),
           nv_state_backward_(N_VMake_Serial(N, state_backward.data())),
           nv_quad_(N_VMake_Serial(num_args_vars, quad.data())),
           nv_absolute_tolerance_forward_(
               N_VMake_Serial(N, absolute_tolerance_forward.data())),
           nv_absolute_tolerance_backward_(
               N_VMake_Serial(N, absolute_tolerance_backward.data())),
           A_forward_(SUNDenseMatrix(N, N)),
           A_backward_(SUNDenseMatrix(N, N)),
           LS_forward_(
               N == 0 ? nullptr
                      : SUNDenseLinearSolver(nv_state_forward_, A_forward_)),
           LS_backward_(
               N == 0 ? nullptr
                      : SUNDenseLinearSolver(nv_state_backward_, A_backward_)),
           N_(N),
           function_name_str_(function_name),
           cvodes_mem_(CVodeCreate(solver_forward)),
           f_(std::forward<FF>(f)) {
       if (cvodes_mem_ == nullptr) {
         throw std::runtime_error("CVodeCreate failed to allocate memory");
       }
     }

     virtual ~cvodes_solver() {
       SUNMatDestroy(A_forward_);
       SUNMatDestroy(A_backward_);
       if (N_ > 0) {
         SUNLinSolFree(LS_forward_);
         SUNLinSolFree(LS_backward_);
       }
       N_VDestroy_Serial(nv_state_forward_);
       N_VDestroy_Serial(nv_state_backward_);
       N_VDestroy_Serial(nv_quad_);
       N_VDestroy_Serial(nv_absolute_tolerance_forward_);
       N_VDestroy_Serial(nv_absolute_tolerance_backward_);

       CVodeFree(&cvodes_mem_);
     }
   };

template <typename ReturnType, typename F, typename T_y0, typename T_t0, typename T_ts,
            typename... T_Args>
  class cvodes_integrator_adjoint_impl<require_arithmetic_t<ReturnType>,
    F, T_y0, T_t0, T_ts, T_Args...> :  {
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

    cvodes_solver<T_Return> solver_{nullptr};
    template <typename ReturnType, typename T_y0_t0, typename T_t0, typename T_ts,
                typename... T_Args>
    struct cvodes_integrator_adjoint_data<void, T_Return,
      T_y0, T_t0, T_ts, T_Args...> cvode_data;
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
    template <typename FF, require_eigen_col_vector_t<T_y0>* = nullptr>
    cvodes_integrator_adjoint_impl(
        const char* function_name, FF&& f, const T_y0& y0, const T_t0& t0,
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
        : cvode_data_(/* Pass members */),
         solver_(function_name, f, N_, num_args_vars_, ts_.size(), solver_forward_,
          state_forward_, state_backward_, quad_, absolute_tolerance_forward_,
          absolute_tolerance_backward_) {

      check_flag_sundials(
          CVodeInit(solver_->cvodes_mem_, &cvodes_integrator_adjoint_impl::cv_rhs,
                    value_of(t0_), solver_->nv_state_forward_),
          "CVodeInit");

      // Assign pointer to this as user data
      check_flag_sundials(
          CVodeSetUserData(solver_->cvodes_mem_, reinterpret_cast<void*>(this)),
          "CVodeSetUserData");

      cvodes_set_options(solver_->cvodes_mem_, max_num_steps_);

      check_flag_sundials(
          CVodeSVtolerances(solver_->cvodes_mem_, relative_tolerance_forward_,
                            solver_->nv_absolute_tolerance_forward_),
          "CVodeSVtolerances");

      check_flag_sundials(
          CVodeSetLinearSolver(solver_->cvodes_mem_, solver_->LS_forward_,
                               solver_->A_forward_),
          "CVodeSetLinearSolver");

      check_flag_sundials(
          CVodeSetJacFn(solver_->cvodes_mem_,
                        &cvodes_integrator_adjoint_impl::cv_jacobian_rhs_states),
          "CVodeSetJacFn");

      /**
       * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
       * times, { t1, t2, t3, ... } using the requested forward solver of CVODES.
       */
      const auto ts_dbl = value_of(ts_);

      double t_init = value_of(t0_);
      for (size_t n = 0; n < ts_dbl.size(); ++n) {
        double t_final = ts_dbl[n];
        if (t_final != t_init) {
            int error_code
                = CVode(solver_->cvodes_mem_, t_final, solver_->nv_state_forward_,
                        &t_init, CV_NORMAL);

            if (unlikely(error_code == CV_TOO_MUCH_WORK)) {
              throw_domain_error(solver_->function_name_str_.c_str(), "", t_final,
                                 "Failed to integrate to next output time (",
                                 ") in less than max_num_steps steps");
            } else {
              check_flag_sundials(error_code, "CVode");
            }
        }
        store_state(n, state_forward_, solver_->y_return_[n]);

        t_init = t_final;
      }
    }

   private:
    void store_state(std::size_t n, const Eigen::VectorXd& state,
                     Eigen::Matrix<double, Eigen::Dynamic, 1>& state_return) {
      y_[n] = state;
      state_return = state;
    }

   public:
    /**
     * Obtain solution of ODE.
     *
     * @return std::vector of Eigen::Matrix of the states of the ODE, one for each
     *   solution time (excluding the initial state)
     */
    std::vector<Eigen::Matrix<T_Return, Eigen::Dynamic, 1>> solution() noexcept {
      return solver_->y_return_;
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
    constexpr static cvodes_integrator_adjoint_impl* cast_to_self(void* mem) {
      return static_cast<cvodes_integrator_adjoint_impl*>(mem);
    }

    /**
     * Calculates the ODE RHS, dy_dt, using the user-supplied functor at
     * the given time t and state y.
     */
    inline int rhs(double t, const double* y, double*& dy_dt) const {
      const Eigen::VectorXd y_vec = Eigen::Map<const Eigen::VectorXd>(y, N_);
      const Eigen::VectorXd dy_dt_vec = rhs(t, y_vec, value_of_args_tuple_);
      check_size_match(solver_->function_name_str_.c_str(), "dy_dt",
                       dy_dt_vec.size(), "states", N_);
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
          = rhs(t, y_var, value_of_args_tuple_);

      check_size_match(solver_->function_name_str_.c_str(), "dy_dt",
                       fy_var.size(), "states", N_);

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
     * @param[in] y State of system
     * @param[in] t Time
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

  /**
   * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
   * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
   * BDF solver or the non-stiff Adams solver from CVODES. The ODE system is
   * integrated using the adjoint sensitivity approach of CVODES. This
   * implementation handles the case of a double return type which ensures that no
   * resources are left on the AD stack.
   *
   * \p f must define an operator() with the signature as:
   *   template<typename T_t, typename T_y, typename... T_Args>
   *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
   *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
   *     std::ostream* msgs, const T_Args&... args);
   *
   * t is the time, y is the state, msgs is a stream for error messages, and args
   * are optional arguments passed to the ODE solve function (which are passed
   * through to \p f without modification).
   *
   * @tparam F Type of ODE right hand side
   * @tparam T_y0 Type of initial state
   * @tparam T_t0 Type of initial time
   * @tparam T_ts Type of output times
   * @tparam T_Args Types of pass-through parameters
   *
   * @param function_name Calling function name (for printing debugging messages)
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
   * @param num_steps_between_checkpoints Number of integrator steps after which a
   * checkpoint is stored for the backward pass
   * @param interpolation_polynomial type of polynomial used for interpolation
   * @param solver_forward solver used for forward pass
   * @param solver_backward solver used for backward pass
   * @param[in, out] msgs the print stream for warning messages
   * @param args Extra arguments passed unmodified through to ODE right hand side
   * @return An `std::vector` of Eigen column vectors with scalars equal to
   *  the least upper bound of `T_y0`, `T_t0`, `T_ts`, and the lambda's arguments.
   *  This represents the solution to ODE at times \p ts
   */
  template <typename F, typename T_y0, typename T_t0, typename T_ts,
            typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
            require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd,
                                           T_abs_tol_bwd>* = nullptr,
            require_all_st_arithmetic<T_y0, T_t0, T_ts, T_Args...>* = nullptr>
  std::vector<Eigen::VectorXd> ode_adjoint_impl(
      const char* function_name, F&& f, const T_y0& y0, const T_t0& t0,
      const std::vector<T_ts>& ts, double relative_tolerance_forward,
      const T_abs_tol_fwd& absolute_tolerance_forward,
      double relative_tolerance_backward,
      const T_abs_tol_bwd& absolute_tolerance_backward,
      double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
      long int max_num_steps,                  // NOLINT(runtime/int)
      long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
      int interpolation_polynomial, int solver_forward, int solver_backward,
      std::ostream* msgs, const T_Args&... args) {
      using T_Return = return_type_t<T_y0, T_t0, T_ts, T_Args...>;

      using cvode_integrator
          = cvodes_integrator_adjoint_impl<void, T_Return, F, plain_type_t<T_y0>, T_t0, T_ts,
                                           plain_type_t<T_Args>...>;

      return cvode_integrator(
          function_name, std::forward<F>(f), y0, t0, ts,
          relative_tolerance_forward, absolute_tolerance_forward,
          relative_tolerance_backward, absolute_tolerance_backward,
          relative_tolerance_quadrature, absolute_tolerance_quadrature,
          max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
          solver_forward, solver_backward, msgs, args...).solution();
  }

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver or the non-stiff Adams solver from CVODES. The ODE system is
 * integrated using the adjoint sensitivity approach of CVODES.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the state, msgs is a stream for error messages, and args
 * are optional arguments passed to the ODE solve function (which are passed
 * through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
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
 * @param num_steps_between_checkpoints Number of integrator steps after which a
 * checkpoint is stored for the backward pass
 * @param interpolation_polynomial type of polynomial used for interpolation
 * @param solver_forward solver used for forward pass
 * @param solver_backward solver used for backward pass
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return An `std::vector` of Eigen column vectors with scalars equal to
 *  the least upper bound of `T_y0`, `T_t0`, `T_ts`, and the lambda's arguments.
 *  This represents the solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd,
                                         T_abs_tol_bwd>* = nullptr>
auto ode_adjoint_tol_ctl(
    F&& f, const T_y0& y0, const T_t0& t0, const std::vector<T_ts>& ts,
    double relative_tolerance_forward,
    const T_abs_tol_fwd& absolute_tolerance_forward,
    double relative_tolerance_backward,
    const T_abs_tol_bwd& absolute_tolerance_backward,
    double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
    long int max_num_steps,                  // NOLINT(runtime/int)
    long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
    int interpolation_polynomial, int solver_forward, int solver_backward,
    std::ostream* msgs, const T_Args&... args) {
  return ode_adjoint_impl(
      "ode_adjoint_tol_ctl", std::forward<F>(f), y0, t0, ts,
      relative_tolerance_forward, absolute_tolerance_forward,
      relative_tolerance_backward, absolute_tolerance_backward,
      relative_tolerance_quadrature, absolute_tolerance_quadrature,
      max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
      solver_forward, solver_backward, msgs, args...);
}

}  // namespace math
}  // namespace stan
#endif
