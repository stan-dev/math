#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_fho.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

using forced_harm_osc_ts_test_types = ::testing::Types<
  std::tuple<ode_rk45_functor, ode_rk45_functor, stan::math::var, double, double>,
  std::tuple<ode_ckrk_functor, ode_ckrk_functor, stan::math::var, double, double>,
  std::tuple<ode_bdf_functor, ode_bdf_functor, stan::math::var, double, double>,
  std::tuple<ode_adams_functor, ode_adams_functor, stan::math::var, double, double>,
  std::tuple<ode_rk45_functor, ode_rk45_functor, stan::math::var, stan::math::var, double>,
  std::tuple<ode_ckrk_functor, ode_ckrk_functor, stan::math::var, stan::math::var, double>,
  std::tuple<ode_bdf_functor, ode_bdf_functor, stan::math::var, stan::math::var, double>,
  std::tuple<ode_adams_functor, ode_adams_functor, stan::math::var, stan::math::var, double>,
  std::tuple<ode_rk45_functor, ode_rk45_functor, stan::math::var, stan::math::var, stan::math::var>,
  std::tuple<ode_ckrk_functor, ode_ckrk_functor, stan::math::var, stan::math::var, stan::math::var>,
  std::tuple<ode_bdf_functor, ode_bdf_functor, stan::math::var, stan::math::var, stan::math::var>,
  std::tuple<ode_adams_functor, ode_adams_functor, stan::math::var, stan::math::var, stan::math::var>>;

TYPED_TEST_SUITE_P(forced_harm_osc_ts_test);
TYPED_TEST_P(forced_harm_osc_ts_test, ts_ad) {
  this->test_ts_ad();
}
REGISTER_TYPED_TEST_SUITE_P(forced_harm_osc_ts_test, ts_ad);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, forced_harm_osc_ts_test,
                               forced_harm_osc_ts_test_types);
