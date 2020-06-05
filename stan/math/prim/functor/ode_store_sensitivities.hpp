#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Build output for state of ODE solve from the current coupled_state when
 * all arguments are arithmetic.
 *
 * @tparam F Type of ODE right hand side
 *
 * @param f Right hand side of the ODE
 * @param coupled_state Current state of the coupled_ode_system
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return ODE state
 */
template <typename F, typename T_y0_t0, typename T_t0, typename T_t,
          typename... Args,
          typename
          = require_all_arithmetic_t<typename F::captured_scalar_t__, T_y0_t0,
                                     T_t0, T_t, scalar_type_t<Args>...>>
Eigen::VectorXd ode_store_sensitivities(
    const F& f, const Eigen::VectorXd& coupled_state,
    const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1>& y0, T_t0 t0, T_t t,
    std::ostream* msgs, const Args&... args) {
  return coupled_state.head(y0.size());
}

}  // namespace math
}  // namespace stan
#endif
