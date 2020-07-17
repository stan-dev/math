#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * When all arguments are arithmetic, there are no sensitivities to store, so
 *  the function just returns the current coupled_state.
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
 * @return ODE state
 */
template <
    typename F, typename T_y0_t0, typename T_t0, typename T_t, typename... Args,
    typename
    = require_all_arithmetic_t<T_y0_t0, T_t0, T_t, scalar_type_t<Args>...>>
Eigen::VectorXd ode_store_sensitivities(
    const F& f, const std::vector<double>& coupled_state,
    const Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1>& y0, T_t0 t0, T_t t,
    std::ostream* msgs, const Args&... args) {
  return Eigen::Map<const Eigen::VectorXd>(coupled_state.data(),
                                           coupled_state.size());
}

}  // namespace math
}  // namespace stan
#endif
