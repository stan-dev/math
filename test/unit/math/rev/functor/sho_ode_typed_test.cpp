#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
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
using harmonic_oscillator_fd_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<ode_adams_functor, ode_bdf_functor, ode_ckrk_functor,
                     ode_rk45_functor, ode_adjoint_functor>,
    ::testing::Types<double>,  // t
    ::testing::Types<double>,  // y0
    ::testing::Types<double>   // theta
    >;

TYPED_TEST_SUITE_P(harmonic_oscillator_test);
TYPED_TEST_P(harmonic_oscillator_test, param_and_data_finite_diff) {
  this->t0 = 0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = 1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = -1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_test,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_test,
                               harmonic_oscillator_fd_test_types);

TYPED_TEST_SUITE_P(harmonic_oscillator_data_test);
TYPED_TEST_P(harmonic_oscillator_data_test, param_and_data_finite_diff) {
  this->t0 = 0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = 1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);

  this->t0 = -1.0;
  for (size_t i = 0; i < this->ts.size(); ++i) {
    this->ts[i] = this->t0 + 0.1 * (i + 1);
  }
  this->test_fd_vd(1.e-8, 1e-4);
  this->test_fd_dv(1.e-8, 1e-4);
  this->test_fd_vv(1.e-8, 1e-4);
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_data_test,
                            param_and_data_finite_diff);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_data_test,
                               harmonic_oscillator_fd_test_types);

using harmonic_oscillator_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<ode_adams_functor, ode_bdf_functor, ode_ckrk_functor,
                     ode_rk45_functor, ode_adjoint_functor>,
    ::testing::Types<double>,                                  // t
    ::testing::Types<double, stan::math::var_value<double> >,  // y0
    ::testing::Types<double, stan::math::var_value<double> >   // theta
    >;

TYPED_TEST_SUITE_P(harmonic_oscillator_t0_ad_test);
TYPED_TEST_P(harmonic_oscillator_t0_ad_test, t0_ad) {
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_rk45_functor>::value) {
    this->test_t0_ad(5e-6);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_ckrk_functor>::value) {
    this->test_t0_ad(5e-6);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_adams_functor>::value) {
    this->test_t0_ad(1e-8);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_bdf_functor>::value) {
    this->test_t0_ad(1e-7);
  }
  if (std::is_same<std::tuple_element_t<0, TypeParam>,
                   ode_adjoint_functor>::value) {
    this->test_t0_ad(1e-7);
  }
}
REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_t0_ad_test, t0_ad);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_t0_ad_test,
                               harmonic_oscillator_test_types);

TYPED_TEST_SUITE_P(harmonic_oscillator_analytical_test);
TYPED_TEST_P(harmonic_oscillator_analytical_test, dv) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_dv.ts = std::vector<double>{t};
  this->ode_sol_dv.y0[0] = chi;
  Eigen::VectorXd x(1);
  x << omega;
  this->test_analytical(this->ode_sol_dv, this->analy_sol_functor(),
                        this->analy_grad_omega_sol_functor(), 1e-5, x, t, omega,
                        chi);
}
TYPED_TEST_P(harmonic_oscillator_analytical_test, vd) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_vd.ts = std::vector<double>{t};
  this->ode_sol_vd.theta[0] = omega;
  Eigen::VectorXd x(1);
  x << chi;
  this->test_analytical(this->ode_sol_vd, this->analy_sol_functor(),
                        this->analy_grad_chi_sol_functor(), 1e-5, x, t, omega,
                        chi);
}

TYPED_TEST_P(harmonic_oscillator_analytical_test, vv) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_vv.ts = std::vector<double>{t};
  Eigen::VectorXd x(2);
  x << omega, chi;
  this->test_analytical(this->ode_sol_vv, this->analy_sol_functor(),
                        this->analy_grad_sol_functor(), 1e-5, x, t, omega, chi);
}

REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_analytical_test, dv, vd, vv);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_analytical_test,
                               harmonic_oscillator_test_types);
