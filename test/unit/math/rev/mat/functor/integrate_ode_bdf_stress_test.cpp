#include <stan/math/rev/mat.hpp>
#include <boost/random.hpp>
#include <stan/math/prim/scal/prob/lognormal_rng.hpp>
#include <gtest/gtest.h>
// very small michaelis menten example
#include <test/unit/math/rev/arr/functor/coupled_mm.hpp>
#include <test/unit/util.hpp>
#include <vector>

// test which triggers the too much work exception from odeint
TEST(StanOde_tooMuchWork_test, cvodes_coupled_mm) {
  coupled_mm_ode_fun f_;

  boost::ecuyer1988 rng;

  // initial value and parameters from model definition

  double t0 = 0;

  std::vector<double> ts_long;
  ts_long.push_back(1E10);

  std::vector<double> ts_short;
  ts_short.push_back(1);

  std::vector<double> data;

  std::vector<int> data_int;

  for (std::size_t i = 0; i < 1000; i++) {
    stan::math::start_nested();

    std::vector<stan::math::var> theta_v(4);

    theta_v[0] = stan::math::lognormal_rng(1, 2, rng);
    theta_v[1] = stan::math::lognormal_rng(-1, 2, rng);
    theta_v[2] = stan::math::lognormal_rng(-1, 2, rng);
    theta_v[3] = stan::math::lognormal_rng(-2, 2, rng);

    std::vector<stan::math::var> y0_v(2);
    y0_v[0] = stan::math::lognormal_rng(5, 2, rng);
    y0_v[1] = stan::math::lognormal_rng(-1, 2, rng);

    std::vector<std::vector<stan::math::var> > res
        = stan::math::integrate_ode_bdf(f_, y0_v, t0, ts_long, theta_v, data,
                                        data_int, 0, 1E-10, 1E-10, 10000);

    stan::math::grad(res[0][0].vi_);

    stan::math::recover_memory_nested();
  }
}
