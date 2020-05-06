#ifndef STAN_MATH_REV_FUNCTOR_ODE_ADD_TIME_GRADIENTS_HPP
#define STAN_MATH_REV_FUNCTOR_ODE_ADD_TIME_GRADIENTS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {
namespace internal {
inline double ode_time_gradient_var(double t, double grad) {
  invalid_argument("ode_time_gradient_var", "t", t,
                   "t must be of type var, not double");
  return 0.0;
}

inline var ode_time_gradient_var(var t, double grad) {
  return var(new precomp_v_vari(0.0, t.vi_, grad));
}
};  // namespace internal

template <typename F, typename T_initial, typename T_out, typename... T_Args>
inline auto ode_add_time_gradients(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
                            double t0, const std::vector<double>& ts,
                            const std::vector<Eigen::Matrix<T_out, Eigen::Dynamic, 1>>& y,
                            std::ostream* msgs, const T_Args&... args) {
  return y;
}

template <typename F, typename T_initial, typename T_t0, typename T_ts,
          typename T_out, typename... T_Args>
std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>>
ode_add_time_gradients(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
		       const T_t0& t0, const std::vector<T_ts>& ts,
		       const std::vector<Eigen::Matrix<T_out, Eigen::Dynamic, 1>>& y,
		       std::ostream* msgs, const T_Args&... args) {
  Eigen::VectorXd f_y0_t0;
  if (is_var<T_t0>::value) {
    f_y0_t0 = f(value_of(t0), value_of(y0), msgs, value_of(args)...);

    check_size_match("coupled_ode_observer", "dy_dt", f_y0_t0.size(), "states",
                     y0.size());
  }

  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y_with_gradients(y.size());

  for (size_t i = 0; i < y.size(); ++i) {
    Eigen::VectorXd f_y_t;

    y_with_gradients[i].resize(y[i].size());

    if (is_var<T_ts>::value)
      f_y_t = f(value_of(ts[i]), value_of(y[i]), msgs, value_of(args)...);

    for (size_t j = 0; j < y[i].size(); ++j) {
      y_with_gradients[i](j) = y[i][j];

      if (is_var<T_t0>::value)
        y_with_gradients[i](j) += internal::ode_time_gradient_var(t0, -f_y0_t0[j]);

      if (is_var<T_ts>::value)
        y_with_gradients[i](j) += internal::ode_time_gradient_var(ts[i], f_y_t[j]);
    }
  }

  return y_with_gradients;
}

}  // namespace math
}  // namespace stan
#endif
