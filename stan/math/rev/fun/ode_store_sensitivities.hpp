#ifndef STAN_MATH_REV_FUN_ODE_STORE_SENSITIVITIES_HPP
#define STAN_MATH_REV_FUN_ODE_STORE_SENSITIVITIES_HPP

#include <stan/math/prim/fun/ode_store_sensitivities.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename T_initial, typename... Args>
Eigen::Matrix<var, Eigen::Dynamic, 1> ode_store_sensitivities(const Eigen::VectorXd& coupled_state,
							      const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
							      const Args&... args) {
  const size_t N = y0.size();
  const size_t y0_vars = count_vars(y0);
  const size_t args_vars = count_vars(args...);
  Eigen::Matrix<var, Eigen::Dynamic, 1> yt(N);

  for (size_t j = 0; j < N; j++) {
    const size_t total_vars = y0_vars + args_vars;

    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(total_vars);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(total_vars);

    vari** varis_ptr = varis;
    double* partials_ptr = partials;

    // iterate over parameters for each equation
    varis_ptr = save_varis(varis_ptr, y0);
    for (std::size_t k = 0; k < y0_vars; k++) {
      //*varis_ptr = y0[k].vi_;
      *partials_ptr = coupled_state(N + y0_vars * k + j);
      partials_ptr++;
    }

    varis_ptr = save_varis(varis_ptr, args...);
    for (std::size_t k = 0; k < args_vars; k++) {
      // dy[j]_dtheta[k]
      // theta[k].vi_
      *partials_ptr = coupled_state(N + N * y0_vars + N * k + j);
      partials_ptr++;
    }

    yt(j) = new precomputed_gradients_vari(coupled_state(j), total_vars,
					   varis, partials);
  }

  return yt;
}

}  // namespace math
}  // namespace stan

#endif
