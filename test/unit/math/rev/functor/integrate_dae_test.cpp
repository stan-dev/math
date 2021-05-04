#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/functor/integrate_dae.hpp>

#include <nvector/nvector_serial.h>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

struct chemical_kinetics {
  template <typename T0, typename TYY, typename TYP, typename TPAR>
  inline std::vector<stan::return_type_t<TYY, TYP, TPAR>> operator()(
      const T0& t_in, const std::vector<TYY>& yy, const std::vector<TYP>& yp,
      const std::vector<TPAR>& theta, const std::vector<double>& x_r,
      const std::vector<int>& x_i, std::ostream* msgs) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<stan::return_type_t<TYY, TYP, TPAR>> res(3);

    auto yy1 = yy.at(0);
    auto yy2 = yy.at(1);
    auto yy3 = yy.at(2);

    auto yp1 = yp.at(0);
    auto yp2 = yp.at(1);

    auto p1 = theta.at(0);
    auto p2 = theta.at(1);
    auto p3 = theta.at(2);

    res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
    res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
    res[2] = yy1 + yy2 + yy3 - 1.0;

    return res;
  }
};

struct StanIntegrateDAETest : public ::testing::Test {
  chemical_kinetics f;
  std::vector<double> yy0;
  std::vector<double> yp0;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  const std::vector<int> eq_id;
  const double t0;
  std::vector<double> ts;

  void SetUp() { stan::math::recover_memory(); }

  StanIntegrateDAETest()
      : yy0{1.0, 0.0, 0.0},
        yp0{-0.04, 0.04, 0.0},
        theta{0.040, 1.0e4, 3.0e7},
        msgs{0},
        eq_id{1, 1, 0},
        t0(0.0) {
    const size_t nout{4};
    const double h{0.4};
    for (size_t i = 0; i < nout; ++i)
      ts.push_back(h * std::pow(10, i));
  }
};

TEST_F(StanIntegrateDAETest, idas_ivp_system_yy0) {
  using stan::math::idas_forward_system;
  using stan::math::integrate_dae;
  std::vector<std::vector<double>> yy
      = integrate_dae(f, yy0, yp0, t0, ts, theta, x_r, x_i, 1e-4, 1e-8);
  EXPECT_NEAR(0.985172, yy[0][0], 1e-6);
  EXPECT_NEAR(0.0147939, yy[0][2], 1e-6);
  EXPECT_NEAR(0.905521, yy[1][0], 1e-6);
  EXPECT_NEAR(0.0944571, yy[1][2], 1e-6);
}

TEST_F(StanIntegrateDAETest, forward_sensitivity_theta) {
  using stan::math::idas_forward_system;
  using stan::math::idas_integrator;
  using stan::math::integrate_dae;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = to_var(theta);

  std::vector<std::vector<var>> yy
      = integrate_dae(f, yy0, yp0, t0, ts, theta_var, x_r, x_i, 1e-5, 1e-12);
  EXPECT_NEAR(0.985172, value_of(yy[0][0]), 1e-6);
  EXPECT_NEAR(0.0147939, value_of(yy[0][2]), 1e-6);
  EXPECT_NEAR(0.905519, value_of(yy[1][0]), 1e-6);
  EXPECT_NEAR(0.0944588, value_of(yy[1][2]), 1e-6);

  // test derivatives against central difference results
  std::vector<double> g;
  const double h = 1.e-2;
  const std::vector<double> theta1{theta[0] - theta[0] * h, theta[1], theta[2]};
  const std::vector<double> theta2{theta[0] + theta[0] * h, theta[1], theta[2]};
  stan::math::idas_integrator solver(1e-5, 1e-12, 1000);
  idas_forward_system<chemical_kinetics, double, double, double> dae1(
      f, eq_id, yy0, yp0, theta1, x_r, x_i, msgs),
      dae2(f, eq_id, yy0, yp0, theta2, x_r, x_i, msgs);
  std::vector<std::vector<double>> yy1 = solver.integrate(dae1, t0, ts);
  std::vector<std::vector<double>> yy2 = solver.integrate(dae2, t0, ts);

  double yys_finite_diff = (yy2[3][1] - yy1[3][1]) / (2.0 * theta[0] * h);
  stan::math::set_zero_all_adjoints();
  yy[3][1].grad(theta_var, g);
  EXPECT_NEAR(yys_finite_diff, g[0], 1e-6);
}

TEST_F(StanIntegrateDAETest, inconsistent_ic_error) {
  using stan::math::integrate_dae;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  std::vector<var> theta_var = to_var(theta);

  yy0.back() = -0.1;
  EXPECT_THROW_MSG(
      integrate_dae(f, yy0, yp0, t0, ts, theta_var, x_r, x_i, 1e-5, 1e-12),
      std::domain_error, "DAE residual at t0");
}
