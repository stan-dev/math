#ifndef STAN_MATH_REV_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_REV_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Build output vars for a state of the ODE solve, storing the sensitivities
 * precomputed using the forward sensitivity problem in precomputed varis.
 *
 * Any combination of y0, t0, ts, and any of the args arguments can have the
 * var scalar type.
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param coupled_state Current state of the coupled_ode_system
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return ODE state with scalar type var
 */
template <
    typename F, typename T_y0_t0, typename T_t0, typename T_t, typename... Args,
    typename = require_any_autodiff_t<typename F::captured_scalar_t__, T_y0_t0,
                                      T_t0, T_t, scalar_type_t<Args>...>>
Eigen::Matrix<var, Eigen::Dynamic, 1> ode_store_sensitivities(
    const F& f, const Eigen::VectorXd& coupled_state,
    const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1>& y0, const T_t0& t0,
    const T_t& t, std::ostream* msgs, const Args&... args) {
  using F_dbl = typename F::ValueOf__;
  const size_t N = y0.size();
  const size_t y0_vars = count_vars(y0);
  const size_t args_vars = count_vars(args...);
  const size_t t0_vars = count_vars(t0);
  const size_t t_vars = count_vars(t);

  Eigen::Matrix<var, Eigen::Dynamic, 1> yt(N);

  Eigen::VectorXd y = coupled_state.head(N);

  Eigen::VectorXd f_y_t;
  Eigen::VectorXd f_y0_t0;

  if (is_var<T_t>::value || is_var<T_t0>::value) {
    F_dbl f_dbl = f;

    if (is_var<T_t>::value)
      f_y_t = f_dbl(value_of(t), y, msgs, value_of(args)...);

    if (is_var<T_t0>::value)
      f_y0_t0 = f_dbl(value_of(t0), value_of(y0), msgs, value_of(args)...);
  }

  for (size_t j = 0; j < N; j++) {
    const size_t total_vars
        = y0_vars + args_vars + t0_vars + t_vars + f.num_vars__;

    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(total_vars);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(total_vars);

    vari** varis_ptr = varis;
    double* partials_ptr = partials;

    // iterate over parameters for each equation
    varis_ptr = save_varis(varis_ptr, y0);
    for (std::size_t k = 0; k < y0_vars; ++k) {
      //*varis_ptr = y0[k].vi_;
      *partials_ptr = coupled_state(N + y0_vars * k + j);
      partials_ptr++;
    }

    varis_ptr = save_varis(varis_ptr, args...);
    for (std::size_t k = 0; k < args_vars; ++k) {
      // dy[j]_dtheta[k]
      // theta[k].vi_
      *partials_ptr = coupled_state(N + N * y0_vars + N * k + j);
      partials_ptr++;
    }

    f.save_varis(varis_ptr);
    varis_ptr += f.num_vars__;
    for (std::size_t k = 0; k < f.num_vars__; ++k) {
      *partials_ptr
          = coupled_state(N + N * y0_vars + N * args_vars + N * k + j);
      partials_ptr++;
    }

    varis_ptr = save_varis(varis_ptr, t0);
    if (t0_vars > 0) {
      double dyt_dt0 = 0.0;
      for (std::size_t k = 0; k < y0_vars; ++k) {
        // dy[j]_dtheta[k]
        // theta[k].vi_
        dyt_dt0 += -f_y0_t0[k] * coupled_state(N + y0_vars * k + j);
      }
      *partials_ptr = dyt_dt0;
      partials_ptr++;
    }

    varis_ptr = save_varis(varis_ptr, t);
    if (t_vars > 0) {
      *partials_ptr = f_y_t[j];
      partials_ptr++;
    }

    yt(j) = new precomputed_gradients_vari(y(j), total_vars, varis, partials);
  }

  return yt;
}

}  // namespace math
}  // namespace stan

#endif
