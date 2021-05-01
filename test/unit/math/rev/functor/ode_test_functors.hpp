#ifndef STAN_MATH_TEST_ODE_TEST_FUNCTORS_HPP
#define STAN_MATH_TEST_ODE_TEST_FUNCTORS_HPP

#include <stan/math/prim/functor/ode_ckrk.hpp>
#include <stan/math/prim/functor/ode_rk45.hpp>
#include <stan/math/rev/functor/ode_bdf.hpp>
#include <stan/math/rev/functor/ode_adams.hpp>
#include <stan/math/rev/functor/ode_adjoint.hpp>
#include <stan/math/prim/functor/integrate_ode_rk45.hpp>

#define STAN_DEF_ODE_SOLVER_FUNCTOR(solver_name, solver_func)                  \
  struct solver_name##_functor {                                               \
    const std::string functor_name = #solver_name;                             \
                                                                               \
    template <typename F, typename T_y0, typename T_t0, typename T_ts,         \
              typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr> \
    std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,  \
                              Eigen::Dynamic, 1>>                              \
    operator()(const F& f, const T_y0& y0, const T_t0& t0,                     \
               const std::vector<T_ts>& ts, std::ostream* msgs,                \
               const Args&... args) {                                          \
      return solver_func(f, y0, t0, ts, msgs, args...);                        \
    }                                                                          \
                                                                               \
    template <typename F, typename T_y0, typename T_t0, typename T_ts,         \
              typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr> \
    std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,  \
                              Eigen::Dynamic, 1>>                              \
    operator()(const F& f, const T_y0& y0_arg, const T_t0& t0,                 \
               const std::vector<T_ts>& ts, double rtol, double atol,          \
               size_t max_num_steps, std::ostream* msgs,                       \
               const Args&... args) {                                          \
      return solver_func##_tol(f, y0_arg, t0, ts, rtol, atol, max_num_steps,   \
                               msgs, args...);                                 \
    }                                                                          \
  };

#define STAN_DEF_STD_ODE_SOLVER_FUNCTOR(solver_name, solver_func)              \
  struct solver_name##_functor {                                               \
    template <typename F, typename T_y0, typename T_param, typename T_t0,      \
              typename T_ts>                                                   \
    std::vector<std::vector<stan::return_type_t<T_y0, T_param, T_t0, T_ts>>>   \
    operator()(const F& f, const std::vector<T_y0>& y0, const T_t0& t0,        \
               const std::vector<T_ts>& ts, const std::vector<T_param>& theta, \
               const std::vector<double>& x, const std::vector<int>& x_int,    \
               std::ostream* msgs = nullptr, double rtol = 1e-10,              \
               double atol = 1e-10, size_t max_num_step = 1e8) {               \
      return solver_func(f, y0, t0, ts, theta, x, x_int, msgs, rtol, atol,     \
                         max_num_step);                                        \
    }                                                                          \
  };

STAN_DEF_ODE_SOLVER_FUNCTOR(ode_adams, stan::math::ode_adams);
STAN_DEF_ODE_SOLVER_FUNCTOR(ode_ckrk, stan::math::ode_ckrk);
STAN_DEF_ODE_SOLVER_FUNCTOR(ode_bdf, stan::math::ode_bdf);
STAN_DEF_ODE_SOLVER_FUNCTOR(ode_rk45, stan::math::ode_rk45);

STAN_DEF_STD_ODE_SOLVER_FUNCTOR(integrate_ode_adams,
                                stan::math::integrate_ode_adams);
STAN_DEF_STD_ODE_SOLVER_FUNCTOR(integrate_ode_bdf,
                                stan::math::integrate_ode_bdf);
STAN_DEF_STD_ODE_SOLVER_FUNCTOR(integrate_ode_rk45,
                                stan::math::integrate_ode_rk45);

struct ode_adjoint_functor {
  const std::string functor_name = "ode_adjoint";

  template <typename F, typename T_y0, typename T_t0, typename T_ts,
            typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr>
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                            Eigen::Dynamic, 1>>
  operator()(const F& f, const T_y0& y0, const T_t0& t0,
             const std::vector<T_ts>& ts, std::ostream* msgs,
             const Args&... args) {
    return (*this)(f, y0, t0, ts, 1E-10, 1E-10, 1000000, msgs, args...);
  }

  template <typename F, typename T_y0, typename T_t0, typename T_ts,
            typename... Args, stan::require_eigen_vector_t<T_y0>* = nullptr>
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                            Eigen::Dynamic, 1>>
  operator()(const F& f, const T_y0& y0_arg, const T_t0& t0,
             const std::vector<T_ts>& ts, double relative_tolerance,
             double absolute_tolerance, size_t max_num_steps,
             std::ostream* msgs, const Args&... args) {
    const int N = y0_arg.size();
    const double relative_tolerance_forward = relative_tolerance / 8.0;
    const double relative_tolerance_backward = relative_tolerance / 4.0;
    const double relative_tolerance_quadrature = relative_tolerance;
    const Eigen::VectorXd absolute_tolerance_forward
        = Eigen::VectorXd::Constant(N, absolute_tolerance / 6.0);
    const Eigen::VectorXd absolute_tolerance_backward
        = Eigen::VectorXd::Constant(N, absolute_tolerance / 3.0);
    const double absolute_tolerance_quadrature = absolute_tolerance;
    const long int num_steps_between_checkpoints = 150;  // NOLINT(runtime/int)
    const int interpolation_polynomial = CV_HERMITE;
    const int solver_forward = CV_BDF;
    const int solver_backward = CV_ADAMS;

    return stan::math::ode_adjoint_tol_ctl(
        f, y0_arg, t0, ts, relative_tolerance_forward,
        absolute_tolerance_forward, relative_tolerance_backward,
        absolute_tolerance_backward, relative_tolerance_quadrature,
        absolute_tolerance_quadrature, max_num_steps,
        num_steps_between_checkpoints, interpolation_polynomial, solver_forward,
        solver_backward, msgs, args...);
  }

  template <typename F, typename T_y0, typename T_t0, typename T_ts,
            typename... T_Args,
            stan::require_eigen_col_vector_t<T_y0>* = nullptr>
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                            Eigen::Dynamic, 1>>
  operator()(const F& f, const T_y0& y0, const T_t0& t0,
             const std::vector<T_ts>& ts, double relative_tolerance_forward,
             const Eigen::VectorXd& absolute_tolerance_forward,
             double relative_tolerance_backward,
             const Eigen::VectorXd& absolute_tolerance_backward,
             double relative_tolerance_quadrature,
             double absolute_tolerance_quadrature,
             long int max_num_steps,                  // NOLINT(runtime/int)
             long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
             int interpolation_polynomial, int solver_forward,
             int solver_backward, std::ostream* msgs, const T_Args&... args) {
    return stan::math::ode_adjoint_tol_ctl(
        f, y0, t0, ts, relative_tolerance_forward, absolute_tolerance_forward,
        relative_tolerance_backward, absolute_tolerance_backward,
        relative_tolerance_quadrature, absolute_tolerance_quadrature,
        max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
        solver_forward, solver_backward, msgs, args...);
  }
};

#endif
