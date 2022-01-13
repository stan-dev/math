#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_dae_index_3.hpp>
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
using dae_test_types = boost::mp11::mp_product<
    ode_test_tuple, ::testing::Types<dae_functor>,
    ::testing::Types<double>,                                  // t
    ::testing::Types<double, stan::math::var_value<double> >,  // yy0
    ::testing::Types<double, stan::math::var_value<double> >,  // yp0
    ::testing::Types<double, stan::math::var_value<double> >   // theta
    >;

TYPED_TEST_SUITE_P(index_3_dae_test);
TYPED_TEST_P(index_3_dae_test, solver_failure) {
  EXPECT_THROW_MSG(this->apply_solver(), std::runtime_error,
                   "Error test failures occurred too many times");
  this->ts = {0.0001};
  EXPECT_THROW_MSG(this->apply_solver_tol(), std::runtime_error,
                   "Error test failures occurred too many times");
}

REGISTER_TYPED_TEST_SUITE_P(index_3_dae_test, solver_failure);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, index_3_dae_test, dae_test_types);
