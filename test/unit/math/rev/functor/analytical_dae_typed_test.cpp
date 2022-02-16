#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_dae_analytical.hpp>
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

TYPED_TEST_SUITE_P(analytical_dae_test);
TYPED_TEST_P(analytical_dae_test, dv) {
  double k = 0.5;
  double t = 5;
  this->dae_sol_dv.ts = std::vector<double>{t};
  Eigen::VectorXd x(1);
  x << k;
  this->test_analytical(this->dae_sol_dv, this->analy_sol_functor(),
                        this->analy_grad_theta_functor(), 1e-8, x, t, k);
}

TYPED_TEST_P(analytical_dae_test, vd) {
  double k = 0.5;
  double t = 5;
  this->dae_sol_vd.ts = std::vector<double>{t};
  this->dae_sol_vd.theta = k;
  Eigen::VectorXd x(3);
  x << 1.0, 0.0, 0.0;
  this->dae_sol_vd.yy0 = x;
  this->test_analytical(this->dae_sol_vd, this->analy_sol_functor(),
                        this->analy_grad_yy0_functor(), 1e-8, x, t, k);
}

REGISTER_TYPED_TEST_SUITE_P(analytical_dae_test, dv, vd);
INSTANTIATE_TYPED_TEST_SUITE_P(StanDAE, analytical_dae_test, dae_test_types);
