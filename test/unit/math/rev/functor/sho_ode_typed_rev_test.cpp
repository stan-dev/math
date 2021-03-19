#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

using harmonic_oscillator_test_types = ::testing::Types<
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, double, double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, double, double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, double, double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, double, double>>;

TYPED_TEST_SUITE_P(harmonic_oscillator_test);
TYPED_TEST_P(harmonic_oscillator_test, param_and_data_finite_diff) {
  this->t0 = 0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = 1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = -1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_test,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_test,
                               harmonic_oscillator_test_types);

TYPED_TEST_SUITE_P(harmonic_oscillator_data_test);
TYPED_TEST_P(harmonic_oscillator_data_test, param_and_data_finite_diff) {
  this->t0 = 0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = 1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = -1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_data_test,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_data_test,
                               harmonic_oscillator_test_types);

using harmonic_oscillator_t0_ad_test_types = ::testing::Types<
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, double, double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, double, double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, double, double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, double, double>,
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, stan::math::var,
               double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, stan::math::var,
               double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, stan::math::var,
               double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, stan::math::var,
               double>>;

TYPED_TEST_SUITE_P(harmonic_oscillator_t0_ad_test);
TYPED_TEST_P(harmonic_oscillator_t0_ad_test, t0_ad) {
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_rk45_functor>::value) {
    this->test_t0_ad(5e-6);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_ckrk_functor>::value) {
    this->test_t0_ad(5e-6);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_adams_functor>::value) {
    this->test_t0_ad(1e-8);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_bdf_functor>::value) {
    this->test_t0_ad(1e-7);
  }
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_t0_ad_test, t0_ad);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_t0_ad_test,
                               harmonic_oscillator_t0_ad_test_types);
