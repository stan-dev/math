#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/functor/gradient.hpp>
#include <test/unit/math/rev/functor/sho_grad_test_fixture.hpp>
#include <test/unit/math/prim/functor/ode_test_integrators.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

TEST_F(sho_grad_ode_test, ckrk_double_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = omega;

  double f;
  Eigen::VectorXd grad(1);

  test_functor_double_var_1<ckrk_integrator> func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-5);

  test_functor_double_var_2<ckrk_integrator> func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-5);
}

TEST_F(sho_grad_ode_test, ckrk_var_double) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = chi;

  double f;
  Eigen::VectorXd grad(1);

  test_functor_var_double_1<ckrk_integrator> func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-5);

  test_functor_var_double_2<ckrk_integrator> func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-5);
}

TEST_F(sho_grad_ode_test, ckrk_var_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(2);
  x(0) = omega;
  x(1) = chi;

  double f;
  Eigen::VectorXd grad(2);

  test_functor_var_var_1<ckrk_integrator> func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-5);
  EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-5);

  test_functor_var_var_2<ckrk_integrator> func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-5);
  EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-5);
}
