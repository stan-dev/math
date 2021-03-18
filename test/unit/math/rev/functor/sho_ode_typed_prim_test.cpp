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
    std::tuple<ode_adams_functor, ode_adams_functor, double, double, double> >;

TYPED_TEST_SUITE_P(harmonic_oscillator_test);
TYPED_TEST_P(harmonic_oscillator_test, no_error) { this->test_good(); }
TYPED_TEST_P(harmonic_oscillator_test, error_conditions) { this->test_bad(); }
TYPED_TEST_P(harmonic_oscillator_test, value) {
  this->test_value(0.0);
  this->test_value(1.0);
  this->test_value(-1.0);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_test, no_error,
                            error_conditions, value);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_test,
                               harmonic_oscillator_test_types);

using harmonic_oscillator_integrate_ode_test_types = ::testing::Types<
    std::tuple<integrate_ode_rk45_functor, integrate_ode_rk45_functor, double,
               double, double>,
    std::tuple<integrate_ode_bdf_functor, integrate_ode_bdf_functor, double,
               double, double>,
    std::tuple<integrate_ode_adams_functor, integrate_ode_adams_functor, double,
               double, double> >;

TYPED_TEST_SUITE_P(harmonic_oscillator_bad_ode_test);
TYPED_TEST_P(harmonic_oscillator_bad_ode_test, bad_ode_error) {
  this->test_bad_ode();
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_bad_ode_test, bad_ode_error);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_bad_ode_test,
                               harmonic_oscillator_integrate_ode_test_types);

TYPED_TEST_SUITE_P(harmonic_oscillator_data_test);
TYPED_TEST_P(harmonic_oscillator_data_test, bad_param_and_data) {
  this->test_bad_param_and_data();
}
TYPED_TEST_P(harmonic_oscillator_data_test, value) {
  this->test_value(0.0);
  this->test_value(1.0);
  this->test_value(-1.0);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_data_test, bad_param_and_data,
                            value);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_data_test,
                               harmonic_oscillator_test_types);
