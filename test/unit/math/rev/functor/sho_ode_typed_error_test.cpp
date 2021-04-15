#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

/**
 * Outer product of test types
 */
using harmonic_oscillator_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<ode_adams_functor, ode_bdf_functor, ode_ckrk_functor,
                     ode_rk45_functor, ode_adjoint_functor>,
    ::testing::Types<double, stan::math::var_value<double>>,  // t
    ::testing::Types<double, stan::math::var_value<double>>,  // y0
    ::testing::Types<double, stan::math::var_value<double>>   // theta
    >;

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

using harmonic_oscillator_integrate_ode_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<integrate_ode_rk45_functor, integrate_ode_bdf_functor,
                     integrate_ode_adams_functor>,
    ::testing::Types<double, stan::math::var_value<double>>,  // t
    ::testing::Types<double, stan::math::var_value<double>>,  // y0
    ::testing::Types<double, stan::math::var_value<double>>   // theta
    >;

TYPED_TEST_SUITE_P(harmonic_oscillator_bad_ode_test);
TYPED_TEST_P(harmonic_oscillator_bad_ode_test, bad_ode_error) {
  this->test_bad_ode();
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_bad_ode_test, bad_ode_error);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_bad_ode_test,
                               harmonic_oscillator_integrate_ode_test_types);
