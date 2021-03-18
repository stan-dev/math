#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_lorenz.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

using lorenz_test_types = ::testing::Types<
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, double, double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, double, double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, double, double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, double, double> >;

TYPED_TEST_SUITE_P(lorenz_test);
TYPED_TEST_P(lorenz_test, param_and_data_finite_diff) {
  this->test_fd_vd(1.e-6, 5e-2);
  this->test_fd_dv(1.e-6, 5e-2);
  this->test_fd_vv(1.e-6, 5e-2);
}
REGISTER_TYPED_TEST_SUITE_P(lorenz_test, param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, lorenz_test, lorenz_test_types);
