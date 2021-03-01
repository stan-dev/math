#ifndef STAN_MATH_TEST_ODE_TEST_INTEGRATORS_HPP
#define STAN_MATH_TEST_ODE_TEST_INTEGRATORS_HPP

#include <stan/math/rev.hpp>
#include <stan/math/prim/functor/ode_ckrk.hpp>
#include <iostream>
#include <sstream>
#include <vector>

struct ckrk_integrator {
  template <typename F, typename T_y0, typename T_t0, typename T_ts,
            typename... Args>
  std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                            Eigen::Dynamic, 1>>
  operator()(const F& f, const T_y0& y0, T_t0 t0, const std::vector<T_ts>& ts,
             std::ostream* msgs, const Args&... args) {
    return stan::math::ode_ckrk(f, y0, t0, ts, msgs, args...);
  }
};

#endif
