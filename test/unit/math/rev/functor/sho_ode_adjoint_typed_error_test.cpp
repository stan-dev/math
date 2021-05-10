#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho_ctl.hpp>
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
using harmonic_oscillator_ctl_test_types = boost::mp11::mp_product<
    ode_test_tuple, ::testing::Types<ode_adjoint_functor>,
    ::testing::Types<double, stan::math::var_value<double>>,  // t
    ::testing::Types<double, stan::math::var_value<double>>,  // y0
    ::testing::Types<double, stan::math::var_value<double>>   // theta
    >;

TYPED_TEST_SUITE_P(harmonic_oscillator_ctl_test);
TYPED_TEST_P(harmonic_oscillator_ctl_test, no_error) { this->test_good(); }
TYPED_TEST_P(harmonic_oscillator_ctl_test, error_conditions) {
  this->test_bad();
}
TYPED_TEST_P(harmonic_oscillator_ctl_test, value) {
  stan::math::nested_rev_autodiff nested;

  this->test_value(0.0);
  this->test_value(1.0);
  this->test_value(-1.0);

  if (std::is_same<std::tuple_element_t<2, TypeParam>, double>::value
      && std::is_same<std::tuple_element_t<3, TypeParam>, double>::value
      && std::is_same<std::tuple_element_t<4, TypeParam>, double>::value) {
    EXPECT_EQ(stan::math::nested_size(), 0);
  } else {
    EXPECT_GT(stan::math::nested_size(), 0);
  }
}

REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_ctl_test, no_error,
                            error_conditions, value);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_ctl_test,
                               harmonic_oscillator_ctl_test_types);
