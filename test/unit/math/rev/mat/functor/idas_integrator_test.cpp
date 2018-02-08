#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/mat/functor/idas_forward_system.hpp>
#include <stan/math/rev/mat/functor/idas_integrator.hpp>
#include <stan/math/rev/mat/functor/integrate_dae.hpp>

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
  inline std::vector<typename stan::return_type<TYY, TYP, TPAR>::type>
  operator()(const T0& t_in, const std::vector<TYY>& yy,
             const std::vector<TYP>& yp, const std::vector<TPAR>& theta,
             const std::vector<double>& x_r, const std::vector<int>& x_i,
             std::ostream* msgs) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<typename stan::return_type<TYY, TYP, TPAR>::type> res(3);

    auto yy1 = yy.at(0);
    auto yy2 = yy.at(1);
    auto yy3 = yy.at(2);

    auto yp1 = yp.at(0);
    auto yp2 = yp.at(1);
    // auto yp3 = yp.at(2);

    auto p1 = theta.at(0);
    auto p2 = theta.at(1);
    auto p3 = theta.at(2);

    res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
    res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
    res[2] = yy1 + yy2 + yy3 - 1.0;

    // jacobian respect to yy

    return res;
  }
};

TEST(IDAS_integrator, idas_ivp_system_consistent_ic) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{-0.04, 0.04, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  const std::vector<int> eq_id{1, 1, 0};
  stan::math::idas_forward_system<chemical_kinetics, double, double, double>
      dae{f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  stan::math::idas_integrator solver(1e-4, 1e-8, 1e6);
  const double t0{0.0};
  const double h{0.4};
  std::vector<double> ts;
  const size_t nout{4};
  for (size_t i = 0; i < nout; ++i)
    ts.push_back(h * std::pow(10, i));

  std::vector<std::vector<double> > yy = solver.integrate(dae, t0, ts);
  EXPECT_NEAR(0.985172, yy[0][0], 1e-6);
  EXPECT_NEAR(0.0147939, yy[0][2], 1e-6);
  EXPECT_NEAR(0.905521, yy[1][0], 1e-6);
  EXPECT_NEAR(0.0944571, yy[1][2], 1e-6);
}

TEST(IDAS_integrator, idas_ivp_system_inconsistent_ic) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{0.1, 0.0, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  const std::vector<int> eq_id{1, 1, 0};
  stan::math::idas_forward_system<chemical_kinetics, double, double, double>
      dae{f, eq_id, yy0, yp0, theta, x_r, x_i, msgs};
  stan::math::idas_integrator solver(1e-4, 1e-8, 1e6);
  const double t0{0.0};
  const double h{0.4};
  std::vector<double> ts;
  const size_t nout{4};
  for (size_t i = 0; i < nout; ++i)
    ts.push_back(h * std::pow(10, i));

  std::vector<std::vector<double> > yy = solver.integrate(dae, t0, ts);
  EXPECT_NEAR(0.985172, yy[0][0], 1e-6);
  EXPECT_NEAR(0.0147939, yy[0][2], 1e-6);
  EXPECT_NEAR(0.905521, yy[1][0], 1e-6);
  EXPECT_NEAR(0.0944571, yy[1][2], 1e-6);
}

TEST(IDAS_integrator, idas_forward_sensitivity_consistent_ic) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{-0.04, 0.04, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  const std::vector<int> eq_id{1, 1, 0};
  stan::math::idas_integrator solver(1e-4, 1e-10, 1e6);
  const double t0{0.0};
  const double h{0.4};
  std::vector<double> ts;
  const size_t nout{4};
  for (size_t i = 0; i < nout; ++i)
    ts.push_back(h * std::pow(10, i));

  {
    std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
    stan::math::idas_forward_system<chemical_kinetics, double, double,
                                    stan::math::var>
        dae{f, eq_id, yy0, yp0, theta_var, x_r, x_i, msgs};
    std::vector<std::vector<stan::math::var> > yy
        = solver.integrate(dae, t0, ts);
    EXPECT_NEAR(0.985172, stan::math::value_of(yy[0][0]), 1e-6);
    EXPECT_NEAR(0.0147939, stan::math::value_of(yy[0][2]), 1e-6);
    EXPECT_NEAR(0.905520, stan::math::value_of(yy[1][0]), 1e-6);
    EXPECT_NEAR(0.0944571, stan::math::value_of(yy[1][2]), 1e-6);

    std::vector<double> g;
    yy[0][0].grad(theta_var, g);
    EXPECT_NEAR(-0.355942, g[0], 1e-6);
    EXPECT_NEAR(9.54689e-08 * 1e8, g[1] * 1e8, 1e-6);
    EXPECT_NEAR(-1.58393e-11 * 1e11, g[2] * 1e11, 1e-5);
  }

  {
    std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
    stan::math::idas_forward_system<chemical_kinetics, double, double,
                                    stan::math::var>
        dae{f, eq_id, yy0, yp0, theta_var, x_r, x_i, msgs};
    std::vector<std::vector<stan::math::var> > yy
        = solver.integrate(dae, t0, ts);
    std::vector<double> g;
    yy[1][1].grad(theta_var, g);
    EXPECT_NEAR(0.000179251 * 1e4, g[0] * 1e4, 1e-5);
    EXPECT_NEAR(-5.83021e-10 * 1e10, g[1] * 1e10, 1e-5);
    EXPECT_NEAR(-2.76228e-13 * 1e13, g[2] * 1e13, 1e-5);
  }

  {
    std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
    stan::math::idas_forward_system<chemical_kinetics, double, double,
                                    stan::math::var>
        dae{f, eq_id, yy0, yp0, theta_var, x_r, x_i, msgs};
    std::vector<std::vector<stan::math::var> > yy
        = solver.integrate(dae, t0, ts);
    std::vector<double> g;
    yy[2][2].grad(theta_var, g);
    EXPECT_NEAR(4.24741, g[0], 1e-5);
    EXPECT_NEAR(-1.37296e-05 * 1e5, g[1] * 1e5, 1e-5);
    EXPECT_NEAR(2.28835e-09 * 1e9, g[2] * 1e9, 1e-5);
  }
}

TEST(IDAS_integrator, idas_forward_sensitivity_yy0_par) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{-0.04, 0.04, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  std::ostream* msgs = 0;

  const std::vector<int> eq_id{1, 1, 0};
  stan::math::idas_integrator solver(1e-4, 1e-8, 1e6);
  const double t0{0.0};
  const double h{0.4};
  std::vector<double> ts;
  const size_t nout{4};
  for (size_t i = 0; i < nout; ++i)
    ts.push_back(h * std::pow(10, i));

  {
    std::vector<stan::math::var> yy0_var = stan::math::to_var(yy0);
    stan::math::idas_forward_system<chemical_kinetics, stan::math::var, double,
                                    double>
        dae{f, eq_id, yy0_var, yp0, theta, x_r, x_i, msgs};
    std::vector<std::vector<stan::math::var> > yy
        = solver.integrate(dae, t0, ts);
    EXPECT_NEAR(0.985172, stan::math::value_of(yy[0][0]), 1e-6);
    EXPECT_NEAR(0.0147939, stan::math::value_of(yy[0][2]), 1e-6);
    EXPECT_NEAR(0.905520, stan::math::value_of(yy[1][0]), 1e-6);
    EXPECT_NEAR(0.0944571, stan::math::value_of(yy[1][2]), 1e-6);

    std::vector<double> g;
    yy[0][0].grad(yy0_var, g);
    EXPECT_NEAR(0.859846, g[0], 1e-6);
    EXPECT_NEAR(-8.72006e-05 * 1e5, g[1] * 1e5, 1e-5);
  }
}

TEST(IDAS_integrator, idas_dae_solver_yy0) {
  chemical_kinetics f;
  std::vector<double> yy0{1.0, 0.0, 0.0};
  std::vector<double> yp0{-0.04, 0.04, 0.0};
  std::vector<double> theta{0.040, 1.0e4, 3.0e7};
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);

  const std::vector<int> eq_id{1, 1, 0};
  stan::math::idas_integrator solver(1e-4, 1e-8, 1e6);
  const double t0{0.0};
  const double h{0.4};
  std::vector<double> ts;
  const size_t nout{4};
  for (size_t i = 0; i < nout; ++i)
    ts.push_back(h * std::pow(10, i));

  {
    std::vector<stan::math::var> yy0_var = stan::math::to_var(yy0);
    std::vector<std::vector<stan::math::var> > yy = integrate_dae(
        f, eq_id, yy0_var, yp0, t0, ts, theta, x_r, x_i, 1e-4, 1e-8);
    EXPECT_NEAR(0.985172, stan::math::value_of(yy[0][0]), 1e-6);
    EXPECT_NEAR(0.0147939, stan::math::value_of(yy[0][2]), 1e-6);
    EXPECT_NEAR(0.905520, stan::math::value_of(yy[1][0]), 1e-6);
    EXPECT_NEAR(0.0944571, stan::math::value_of(yy[1][2]), 1e-6);

    std::vector<double> g;
    yy[0][0].grad(yy0_var, g);
    EXPECT_NEAR(0.859846, g[0], 1e-6);
    EXPECT_NEAR(-8.72006e-05 * 1e5, g[1] * 1e5, 1e-5);
  }

  {
    std::vector<stan::math::var> yy0_var = stan::math::to_var(yy0);
    std::vector<stan::math::var> yp0_var = stan::math::to_var(yp0);
    std::vector<std::vector<stan::math::var> > yy = integrate_dae(
        f, eq_id, yy0_var, yp0_var, t0, ts, theta, x_r, x_i, 1e-4, 1e-8);
    EXPECT_NEAR(0.985172, stan::math::value_of(yy[0][0]), 1e-6);
    EXPECT_NEAR(0.0147939, stan::math::value_of(yy[0][2]), 1e-6);
    EXPECT_NEAR(0.905520, stan::math::value_of(yy[1][0]), 1e-6);
    EXPECT_NEAR(0.0944571, stan::math::value_of(yy[1][2]), 1e-6);

    std::vector<double> g;
    yy[0][0].grad(yy0_var, g);
    EXPECT_NEAR(0.859846, g[0], 1e-6);
    EXPECT_NEAR(-8.72006e-05 * 1e5, g[1] * 1e5, 1e-5);
  }

  {
    std::vector<stan::math::var> yy0_var = stan::math::to_var(yy0);
    std::vector<stan::math::var> yp0_var = stan::math::to_var(yp0);
    std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
    std::vector<std::vector<stan::math::var> > yy = integrate_dae(
        f, eq_id, yy0_var, yp0_var, t0, ts, theta_var, x_r, x_i, 1e-4, 1e-8);
    EXPECT_NEAR(0.985172, stan::math::value_of(yy[0][0]), 1e-6);
    EXPECT_NEAR(0.0147939, stan::math::value_of(yy[0][2]), 1e-6);
    EXPECT_NEAR(0.905520, stan::math::value_of(yy[1][0]), 1e-6);
    EXPECT_NEAR(0.0944571, stan::math::value_of(yy[1][2]), 1e-6);

    std::vector<double> g;
    yy[0][0].grad(yy0_var, g);
    EXPECT_NEAR(0.859846, g[0], 1e-6);
    EXPECT_NEAR(-8.7197700e-05 * 1e5, g[1] * 1e5, 1e-5);
  }
}

TEST(IDAS_integrator, constructor_error) {
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
}
