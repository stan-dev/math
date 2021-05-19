#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_fho.hpp>
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
using forced_harm_osc_ts_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<ode_adams_functor, ode_bdf_functor, ode_ckrk_functor,
                     ode_rk45_functor, ode_adjoint_functor>,
    ::testing::Types<stan::math::var_value<double> >,          // time
    ::testing::Types<double, stan::math::var_value<double> >,  // y0
    ::testing::Types<double, stan::math::var_value<double> >   // theta
    >;

TYPED_TEST_SUITE_P(forced_harm_osc_ts_test);
TYPED_TEST_P(forced_harm_osc_ts_test, ts_ad) { this->test_ts_ad(); }
REGISTER_TYPED_TEST_SUITE_P(forced_harm_osc_ts_test, ts_ad);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, forced_harm_osc_ts_test,
                               forced_harm_osc_ts_test_types);
