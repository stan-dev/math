#ifndef STAN_MATH_PRIM_FUN_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_PRIM_FUN_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename F, typename... Args,
          typename = require_all_arithmetic_t<scalar_type_t<Args>...>>
Eigen::VectorXd ode_store_sensitivities(const F& f,
                                        const Eigen::VectorXd& coupled_state,
                                        const Eigen::VectorXd& y0, double t0,
                                        double t, std::ostream* msgs,
                                        const Args&... args) {
  return coupled_state.head(y0.size());
}

}  // namespace math
}  // namespace stan
#endif
