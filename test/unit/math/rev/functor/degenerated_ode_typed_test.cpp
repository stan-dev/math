#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_linear_dae.hpp>
#include <test/unit/math/rev/functor/dae_test_functors.hpp>

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

using degenerated_dae_test_types = boost::mp11::mp_product<
    ode_test_tuple, ::testing::Types<dae_functor>,
    ::testing::Types<double>,  // t
    ::testing::Types<double, stan::math::var_value<double>>,
    ::testing::Types<double, stan::math::var_value<double>>,
    ::testing::Types<double, stan::math::var_value<double>>>;
TYPED_TEST_SUITE_P(degenerated_dae_test);
TYPED_TEST_P(degenerated_dae_test, y0_sens) { this->test_ode_sens_y0(); }
TYPED_TEST_P(degenerated_dae_test, theta_sens) { this->test_ode_sens_theta(); }

REGISTER_TYPED_TEST_SUITE_P(degenerated_dae_test, y0_sens, theta_sens);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, degenerated_dae_test,
                               degenerated_dae_test_types);
