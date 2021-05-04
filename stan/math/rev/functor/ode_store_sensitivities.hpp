#ifndef STAN_MATH_REV_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_REV_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Build output vars for a state of the ODE solve, storing the sensitivities
 * precomputed using the forward sensitivity problem in precomputed varis.
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0_t0 Type of initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param coupled_state Current state of the coupled_ode_system
 * @param y0 Initial state
 * @param t0 Initial time
 * @param t Times at which to solve the ODE at
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return ODE state with scalar type var
 */
template <typename F, typename T_y0_t0, typename T_t0, typename T_t,
          typename... Args,
          require_any_autodiff_t<T_y0_t0, T_t0, T_t,
                                 scalar_type_t<Args>...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> ode_store_sensitivities(
    const F& f, const std::vector<double>& coupled_state,
    const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1>& y0, const T_t0& t0,
    const T_t& t, std::ostream* msgs, const Args&... args) {
  const size_t N = y0.size();
  const size_t num_y0_vars = count_vars(y0);
  const size_t num_args_vars = count_vars(args...);
  const size_t num_t0_vars = count_vars(t0);
  const size_t num_t_vars = count_vars(t);
  Eigen::Matrix<var, Eigen::Dynamic, 1> yt(N);

  Eigen::VectorXd y(N);
  for (size_t n = 0; n < N; ++n) {
    y.coeffRef(n) = coupled_state[n];
  }

  Eigen::VectorXd f_y_t;
  if (is_var<T_t>::value)
    f_y_t = f(value_of(t), y, msgs, eval(value_of(args))...);

  Eigen::VectorXd f_y0_t0;
  if (is_var<T_t0>::value)
    f_y0_t0
        = f(value_of(t0), eval(value_of(y0)), msgs, eval(value_of(args))...);

  const size_t total_vars
      = num_y0_vars + num_args_vars + num_t0_vars + num_t_vars;

  vari** varis
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(total_vars);

  save_varis(varis, y0, args..., t0, t);

  // memory for a column major jacobian
  double* jacobian_mem
      = ChainableStack::instance_->memalloc_.alloc_array<double>(N
                                                                 * total_vars);

  Eigen::Map<Eigen::MatrixXd> jacobian(jacobian_mem, total_vars, N);

  for (size_t j = 0; j < N; ++j) {
    for (size_t k = 0; k < num_y0_vars; ++k) {
      jacobian.coeffRef(k, j) = coupled_state[N + num_y0_vars * k + j];
    }

    for (size_t k = 0; k < num_args_vars; ++k) {
      jacobian.coeffRef(num_y0_vars + k, j)
          = coupled_state[N + N * num_y0_vars + N * k + j];
    }

    if (is_var<T_t0>::value) {
      double dyt_dt0 = 0.0;
      for (size_t k = 0; k < num_y0_vars; ++k) {
        dyt_dt0 -= f_y0_t0.coeffRef(k) * coupled_state[N + num_y0_vars * k + j];
      }
      jacobian.coeffRef(num_y0_vars + num_args_vars, j) = dyt_dt0;
    }

    if (is_var<T_t>::value) {
      jacobian.coeffRef(num_y0_vars + num_args_vars + num_t0_vars, j)
          = f_y_t.coeffRef(j);
    }

    // jacobian_mem + j * total_vars points to the jth column of jacobian
    yt(j) = new precomputed_gradients_vari(y(j), total_vars, varis,
                                           jacobian_mem + j * total_vars);
  }

  return yt;
}

}  // namespace math
}  // namespace stan

#endif
