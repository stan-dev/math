#ifndef STAN_MATH_REV_FUNCTOR_COUPLED_ODE_RHS_HPP
#define STAN_MATH_REV_FUNCTOR_COUPLED_ODE_RHS_HPP

#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename F, typename T_y, typename... Args>
std::vector<double> coupled_ode_rhs(const F& f,
				    double t,
				    const std::vector<T_y>& y,
				    std::ostream* msgs,
				    const Args&... args) const {
  // Run nested autodiff in this scope
  nested_rev_autodiff nested;

  int N = y.size();
  const size_t y_vars = count_vars(y);
  const size_t args_vars = count_vars(args...);
  std::vector<double> dz_dt(N + N * y_vars + N * args_vars);
  const std::vector<var> y_vars(y.begin(), y.end());
  std::tuple<decltype(deep_copy_vars(args))...> local_args_tuple(deep_copy_vars(args)...);
    
  std::vector<var> dy_dt_vars = apply([&](auto&&... args) {
      f_(t, y_vars, msgs_, args...);
    }, local_args_tuple);

  check_size_match("coupled_ode_system", "dy_dt", dy_dt_vars.size(), "states",
		   N);

  for (size_t i = 0; i < N; ++i) {
    dz_dt[i] = dy_dt_vars[i].val();
    dy_dt_vars[i].grad();

    for (size_t j = 0; j < y0_vars_; ++j) {
      // orders derivatives by equation (i.e. if there are 2 eqns
      // (y1, y2) and 2 parameters (a, b), dy_dt will be ordered as:
      // dy1_dt, dy2_dt, dy1_da, dy2_da, dy1_db, dy2_db
      double temp_deriv = 0;
      const size_t offset = N + N * j;
      for (size_t k = 0; k < N; k++) {
	temp_deriv += z[N + N * j + k] * y_vars[k].adj();
      }
	
      dz_dt[N + N * j + i] = temp_deriv;
    }
      
    Eigen::VectorXd args_adjoints = Eigen::VectorXd::Zero(args_vars);
    apply([&](auto&&... args) {
	accumulate_adjoints(args_adjoints.data(), args...);
      }, local_args_tuple);

    for (size_t j = 0; j < args_vars; j++) {
      double temp_deriv = args_adjoints(j);
      for (size_t k = 0; k < N; k++) {
	temp_deriv += z[N + N * y0_vars + N * j + k] * y_vars[k].adj();
      }
	
      dz_dt[N + N * y0_vars + N * j + i] = temp_deriv;
    }
      
    nested.set_zero_all_adjoints();
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
