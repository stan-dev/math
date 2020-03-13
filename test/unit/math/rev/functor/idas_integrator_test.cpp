#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/functor/idas_forward_system.hpp>
#include <stan/math/rev/functor/idas_integrator.hpp>

#include <nvector/nvector_serial.h>

#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

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

struct prey_predator_harvest {
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

    constexpr auto r1 = 1.0;
    constexpr auto r2 = 3.0;
    constexpr auto p = 2.0;
    auto mu = theta.at(0);

    res[0] = yp1 - yy1 * (r1 - yy2);
    res[1] = yp2 - yy2 * (r2 - yy2 / yy1 - yy3);
    res[2] = yy3 * (p * yy2 - 1.0) - mu;

    return res;
  }
};

struct IDASIntegratorTest : public ::testing::Test {
  chemical_kinetics f;
  prey_predator_harvest f2;
  std::vector<double> yy0;
  std::vector<double> yp0;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  std::vector<int> eq_id;
  const double t0;
  std::vector<double> ts;

  void SetUp() { stan::math::recover_memory(); }

  IDASIntegratorTest()
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

TEST_F(IDASIntegratorTest, idas_ivp_system_yy0) {
  using stan::math::idas_forward_system;
  using stan::math::idas_integrator;
  idas_forward_system<chemical_kinetics, double, double, double> dae{
      f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  idas_integrator solver(1e-4, 1e-8, 1e6);

  std::vector<std::vector<double>> yy = solver.integrate(dae, t0, ts);
  EXPECT_NEAR(0.985172, yy[0][0], 1e-6);
  EXPECT_NEAR(0.0147939, yy[0][2], 1e-6);
  EXPECT_NEAR(0.905521, yy[1][0], 1e-6);
  EXPECT_NEAR(0.0944571, yy[1][2], 1e-6);
}

TEST_F(IDASIntegratorTest, forward_sensitivity_theta_chemical) {
  using stan::math::idas_forward_system;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  stan::math::idas_integrator solver(1e-5, 1e-12, 1e6);
  std::vector<var> theta_var = to_var(theta);

  idas_forward_system<chemical_kinetics, double, double, var> dae(
      f, eq_id, yy0, yp0, theta_var, x_r, x_i, msgs);
  std::vector<std::vector<var>> yy = solver.integrate(dae, t0, ts);
  EXPECT_NEAR(0.985172, value_of(yy[0][0]), 1e-6);
  EXPECT_NEAR(0.0147939, value_of(yy[0][2]), 1e-6);
  EXPECT_NEAR(0.905519, value_of(yy[1][0]), 1e-6);
  EXPECT_NEAR(0.0944588, value_of(yy[1][2]), 1e-6);

  // test derivatives against central difference results
  std::vector<double> g;
  const double h = 1.e-2;
  const std::vector<double> theta1{theta[0] - theta[0] * h, theta[1], theta[2]};
  const std::vector<double> theta2{theta[0] + theta[0] * h, theta[1], theta[2]};
  idas_forward_system<chemical_kinetics, double, double, double> dae1(
      f, eq_id, yy0, yp0, theta1, x_r, x_i, msgs),
      dae2(f, eq_id, yy0, yp0, theta2, x_r, x_i, msgs);
  std::vector<std::vector<double>> yy1 = solver.integrate(dae1, t0, ts);
  std::vector<std::vector<double>> yy2 = solver.integrate(dae2, t0, ts);

  double yys_finite_diff;
  for (size_t i = 0; i < yy.size(); ++i) {
    yys_finite_diff = (yy2[i][1] - yy1[i][1]) / (2.0 * theta[0] * h);
    stan::math::set_zero_all_adjoints();
    yy[i][1].grad(theta_var, g);
    EXPECT_NEAR(yys_finite_diff, g[0], 5e-6);
  }
}

TEST_F(IDASIntegratorTest, forward_sensitivity_theta_prey) {
  using stan::math::idas_forward_system;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;

  yy0[0] = 0.5;
  yy0[1] = 1.0;
  yy0[2] = 1.0;
  std::fill(yp0.begin(), yp0.end(), 0.0);
  eq_id[0] = 1;
  eq_id[1] = 1;
  for (size_t i = 0; i < ts.size(); ++i)
    ts[i] = (i + 1) * 0.1;
  theta[0] = 1.0;

  stan::math::idas_integrator solver(1e-4, 1e-8, 1000);
  std::vector<var> theta_var = to_var(theta);

  idas_forward_system<prey_predator_harvest, double, double, var> dae(
      f2, eq_id, yy0, yp0, theta_var, x_r, x_i, msgs);
  std::vector<std::vector<var>> yy = solver.integrate(dae, t0, ts);

  // test derivatives against central difference results
  std::vector<double> g;
  const double h = 1.e-5;
  const std::vector<double> theta1{theta[0] - theta[0] * h};
  const std::vector<double> theta2{theta[0] + theta[0] * h};
  idas_forward_system<prey_predator_harvest, double, double, double> dae1(
      f2, eq_id, yy0, yp0, theta1, x_r, x_i, msgs),
      dae2(f2, eq_id, yy0, yp0, theta2, x_r, x_i, msgs);
  std::vector<std::vector<double>> yy1 = solver.integrate(dae1, t0, ts);
  std::vector<std::vector<double>> yy2 = solver.integrate(dae2, t0, ts);

  double yys_finite_diff;
  for (size_t i = 0; i < yy.size(); ++i) {
    yys_finite_diff = (yy2[i][1] - yy1[i][1]) / (2.0 * theta[0] * h);
    stan::math::set_zero_all_adjoints();
    yy[i][1].grad(theta_var, g);
    EXPECT_NEAR(yys_finite_diff, g[0], 1e-6);
  }
}

TEST_F(IDASIntegratorTest, error_handling) {
  using stan::math::idas_forward_system;
  using stan::math::idas_integrator;
  const double rtol = 1e-4;
  const double atol = 1e-8;
  const int n = 600;

  ASSERT_NO_THROW(idas_integrator(rtol, atol, n));

  EXPECT_THROW_MSG(idas_integrator(-1.0E-4, atol, n), std::invalid_argument,
                   "relative tolerance");

  EXPECT_THROW_MSG(idas_integrator(2.E-3, atol, n), std::invalid_argument,
                   "relative tolerance");

  EXPECT_THROW_MSG(idas_integrator(rtol, -1.E-9, n), std::invalid_argument,
                   "absolute tolerance");

  EXPECT_THROW_MSG(idas_integrator(rtol, atol, -100), std::invalid_argument,
                   "max_num_steps");

  idas_forward_system<chemical_kinetics, double, double, double> dae{
      f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  idas_integrator solver(rtol, atol, n);
  double bad_t0 = std::numeric_limits<double>::infinity();
  std::vector<double> bad_ts{std::numeric_limits<double>::infinity()};
  EXPECT_THROW_MSG(solver.integrate(dae, bad_t0, ts), std::domain_error,
                   "initial time");
  EXPECT_THROW_MSG(solver.integrate(dae, t0, bad_ts), std::domain_error,
                   "times");
  bad_t0 = 0;
  bad_ts[0] = -1;
  EXPECT_THROW_MSG(solver.integrate(dae, bad_t0, bad_ts), std::domain_error,
                   "initial time");

  bad_ts[0] = 0.0;
  bad_ts.push_back(0.0);
  EXPECT_THROW_MSG(solver.integrate(dae, t0, bad_ts), std::domain_error,
                   "times");

  bad_ts.clear();
  EXPECT_THROW_MSG(solver.integrate(dae, t0, bad_ts), std::invalid_argument,
                   "times");
}
