#ifndef STAN_MATH_TEST_DAE_TEST_FUNCTORS_HPP
#define STAN_MATH_TEST_DAE_TEST_FUNCTORS_HPP

#include <stan/math/rev/functor/dae.hpp>

struct dae_functor {
  template <typename F, typename T_yy, typename T_yp, typename... T_Args>
  std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
  operator()(const F& f, const T_yy& yy0, const T_yp& yp0, double t0,
             const std::vector<double>& ts, std::ostream* msgs,
             const T_Args&... args) {
    return stan::math::dae(f, yy0, yp0, t0, ts, msgs, args...);
  }

  template <typename F, typename T_yy, typename T_yp, typename... T_Args>
  std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
  operator()(const F& f, const T_yy& yy0, const T_yp& yp0, double t0,
             const std::vector<double>& ts, double rtol, double atol,
             int64_t max_num_steps, std::ostream* msgs, const T_Args&... args) {
    return stan::math::dae_tol(f, yy0, yp0, t0, ts, rtol, atol, max_num_steps,
                               msgs, args...);
  }
};

#endif
