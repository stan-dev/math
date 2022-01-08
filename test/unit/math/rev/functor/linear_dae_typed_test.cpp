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

/**
 * Outer product of test types
 */
using linear_dae_test_types
    = boost::mp11::mp_product<ode_test_tuple, ::testing::Types<dae_functor>,
                              ::testing::Types<double>,  // t
                              ::testing::Types<double>,  // yy0
                              ::testing::Types<double>,  // yp0
                              ::testing::Types<double>   // theta
                              >;

TYPED_TEST_SUITE_P(linear_dae_test);
TYPED_TEST_P(linear_dae_test, analytical_value) {
  double theta = 3.0;
  this->theta = theta;
  Eigen::VectorXd x(1);
  x << theta;
  const double pi = 3.14159265358979323846;
  const double t = 2.3 * pi;
  this->dae_sol.ts = std::vector<double>{t};
  this->test_analytical(this->dae_sol, this->analy_sol_functor(), 1e-6, x, t);
}

TYPED_TEST_P(linear_dae_test, finite_diff) { this->test_fd_dv(1.e-3, 5e-6); }

REGISTER_TYPED_TEST_SUITE_P(linear_dae_test, analytical_value, finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, linear_dae_test, linear_dae_test_types);
