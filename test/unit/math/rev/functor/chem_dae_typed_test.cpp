#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_dae_chem.hpp>
#include <test/unit/math/rev/functor/dae_test_functors.hpp>

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

/**
 * Outer product of test types
 */
using chemical_kinetics_test_types
    = boost::mp11::mp_product<ode_test_tuple, ::testing::Types<dae_functor>,
                              ::testing::Types<double>,  // t
                              ::testing::Types<double>,  // yy0
                              ::testing::Types<double>,  // yp0
                              ::testing::Types<double>   // theta
                              >;

TYPED_TEST_SUITE_P(chemical_kinetics_test);
TYPED_TEST_P(chemical_kinetics_test, param_and_data_finite_diff) {
  // params that gives extreme gradients
  this->test_fd_dv(1.e-3, 2e-3);

  // params that gives moderate gradients
  this->theta = {0.4, 0.5, 0.5};
  this->test_fd_dv(1.e-3, 5e-6);
}
TYPED_TEST_P(chemical_kinetics_test, value) {
  this->test_value(0.0);
  this->test_value(1.0);
  this->test_value(-1.0);
}

TYPED_TEST_P(chemical_kinetics_test, error) { this->test_bad(); }

REGISTER_TYPED_TEST_SUITE_P(chemical_kinetics_test, value, error,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, chemical_kinetics_test,
                               chemical_kinetics_test_types);

TYPED_TEST_SUITE_P(chemical_kinetics_data_test);
TYPED_TEST_P(chemical_kinetics_data_test, value) {
  this->test_value(0.0);
  this->test_value(1.0);
  this->test_value(-1.0);
}
TYPED_TEST_P(chemical_kinetics_data_test, param_and_data_finite_diff) {
  // params that gives extreme gradients
  this->test_fd_dv(1.e-3, 2e-3);

  // params that gives moderate gradients
  this->theta = {0.4, 0.5, 0.5};
  this->x_r = {0.5, 0.5};
  this->test_fd_dv(1.e-3, 2e-6);
}
REGISTER_TYPED_TEST_SUITE_P(chemical_kinetics_data_test, value,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, chemical_kinetics_data_test,
                               chemical_kinetics_test_types);
