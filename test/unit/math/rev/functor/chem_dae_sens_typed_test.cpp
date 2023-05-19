#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_dae_chem.hpp>
#include <test/unit/math/rev/functor/dae_test_functors.hpp>

namespace chem_dae_sens_typed_test {

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

/**
 * Outer product of test types
 */
using chemical_kinetics_sens_test_types = boost::mp11::mp_product<
    ode_test_tuple, ::testing::Types<dae_functor>,
    ::testing::Types<double>,  // t
    ::testing::Types<double,
                     stan::math::var_value<double> >,  // yy
    ::testing::Types<double,
                     stan::math::var_value<double> >,    // yp
    ::testing::Types<stan::math::var_value<double> > >;  // theta

TYPED_TEST_SUITE_P(chemical_kinetics_test);
TYPED_TEST_P(chemical_kinetics_test, value) { this->test_value(0.0); }
TYPED_TEST_P(chemical_kinetics_test, sens) { this->test_sens(0.0); }
REGISTER_TYPED_TEST_SUITE_P(chemical_kinetics_test, value, sens);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, chemical_kinetics_test,
                               chemical_kinetics_sens_test_types);
}  // namespace chem_dae_sens_typed_test
