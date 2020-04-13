#ifndef STAN_MATH_PRIM_FUN_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_PRIM_FUN_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/meta.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename... Args, typename = require_all_arithmetic_t<scalar_type_t<Args>...>>
std::vector<double> ode_store_sensitivities(const std::vector<double>& coupled_state,
					    const std::vector<double>& y0,
					    const Args&... args) {
  return std::vector<double>(coupled_state.data(), coupled_state.data() + y0.size());
}

}  // namespace math
}  // namespace stan
#endif
