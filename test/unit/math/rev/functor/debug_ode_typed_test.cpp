#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

using harmonic_oscillator_test_types = ::testing::Types<
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, double, double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, double, double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, double, double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, double, double>,
  std::tuple<ode_adjoint_functor, ode_adjoint_functor, double, double, double>,
    std::tuple<ode_rk45_functor, ode_rk45_functor, double,  stan::math::var_value<double, void>,
               double>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double,  stan::math::var_value<double, void>,
               double>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double,  stan::math::var_value<double, void>,
               double>,
    std::tuple<ode_adams_functor, ode_adams_functor, double,  stan::math::var_value<double, void>,
               double>,
    std::tuple<ode_adjoint_functor, ode_adjoint_functor, double,  stan::math::var_value<double, void>,
               double>,
    std::tuple<ode_rk45_functor, ode_rk45_functor, double, double,
                stan::math::var_value<double, void>>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double, double,
                stan::math::var_value<double, void>>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double, double,
                stan::math::var_value<double, void>>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, double,
                stan::math::var_value<double, void>>,
    std::tuple<ode_adjoint_functor, ode_adjoint_functor, double, double,
                stan::math::var_value<double, void>>,
    std::tuple<ode_rk45_functor, ode_rk45_functor, double,  stan::math::var_value<double, void>,
                stan::math::var_value<double, void>>,
    std::tuple<ode_ckrk_functor, ode_ckrk_functor, double,  stan::math::var_value<double, void>,
                stan::math::var_value<double, void>>,
    std::tuple<ode_bdf_functor, ode_bdf_functor, double,  stan::math::var_value<double, void>, stan::math::var_value<double, void>>,
    std::tuple<ode_adams_functor, ode_adams_functor, double, stan::math::var_value<double, void>, stan::math::var_value<double, void>>,
  std::tuple<ode_adjoint_functor, ode_adjoint_functor, double, stan::math::var,
               stan::math::var>>;

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
