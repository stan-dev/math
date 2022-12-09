#include <stan/math.hpp>
#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/functor/dae_system.hpp>
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

namespace idas_system_test {

static sundials::Context sundials_context;

struct chemical_kinetics {
  template <typename T0, typename Tyy, typename Typ, typename Tpar>
  inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
      const T0& t_in, const Eigen::Matrix<Tyy, -1, 1>& yy,
      const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
      const std::vector<Tpar>& theta, const std::vector<double>& x_r,
      const std::vector<int>& x_i) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

    auto yy1 = yy(0);
    auto yy2 = yy(1);
    auto yy3 = yy(2);

    auto yp1 = yp(0);
    auto yp2 = yp(1);
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

struct StanIntegrateDAETest : public ::testing::Test {
  chemical_kinetics f;
  Eigen::VectorXd yy0;
  Eigen::VectorXd yp0;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  const std::vector<int> eq_id;
  const double t0;
  std::vector<double> ts;

  void SetUp() { stan::math::recover_memory(); }

  StanIntegrateDAETest()
      : yy0(3),
        yp0(3),
        theta{0.040, 1.0e4, 3.0e7},
        msgs{0},
        eq_id{1, 1, 0},
        t0(0.0) {
    const size_t nout{4};
    const double h{0.4};
    yy0 << 1.0, 0.0, 0.0;
    yp0 << -0.04, 0.04, 0.0;
    for (size_t i = 0; i < nout; ++i)
      ts.push_back(h * std::pow(10, i));
  }
};

TEST_F(StanIntegrateDAETest, dae_system_test) {
  using stan::math::dae_system;

  Eigen::Matrix<stan::math::var, -1, 1> yy0_var = stan::math::to_var(yy0);
  Eigen::Matrix<stan::math::var, -1, 1> yp0_var = stan::math::to_var(yp0);
  std::vector<stan::math::var> theta_var = stan::math::to_var(theta);
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs = nullptr;

  dae_system<chemical_kinetics, Eigen::VectorXd, Eigen::VectorXd,
             std::vector<double>, std::vector<double>, std::vector<int>>
      dae1(f, yy0, yp0, msgs, theta, x_r, x_i);

  // auxiliary variables required due to static member restriction on
  // addressing: "You can take the address of a static member if (and only if)
  // it has an out-of-class definition"
  // For some reason _not_ doing this still works at O3, but not O0.
  auto d1_yy0 = dae1.is_var_yy0;
  EXPECT_FALSE(d1_yy0);
  auto d1_yp0 = dae1.is_var_yp0;
  EXPECT_FALSE(d1_yp0);
  auto d1_par = dae1.is_var_par;
  EXPECT_FALSE(d1_par);
  auto d1_fwd = dae1.use_fwd_sens;
  EXPECT_FALSE(d1_fwd);

  dae_system<chemical_kinetics, Eigen::Matrix<stan::math::var, -1, 1>,
             Eigen::Matrix<stan::math::var, -1, 1>, std::vector<double>>
      dae2(f, yy0_var, yp0_var, msgs, theta);
  auto d2_yy0 = dae2.is_var_yy0;
  EXPECT_TRUE(d2_yy0);
  auto d2_yp0 = dae2.is_var_yp0;
  EXPECT_TRUE(d2_yp0);
  auto d2_par = dae2.is_var_par;
  EXPECT_FALSE(d2_par);
  auto d2_fwd = dae2.use_fwd_sens;
  EXPECT_TRUE(d2_fwd);

  dae_system<chemical_kinetics, Eigen::Matrix<stan::math::var, -1, 1>,
             Eigen::Matrix<stan::math::var, -1, 1>,
             std::vector<stan::math::var>>
      dae3(f, yy0_var, yp0_var, msgs, theta_var);

  auto d3_yy0 = dae3.is_var_yy0;
  EXPECT_TRUE(d3_yy0);
  auto d3_yp0 = dae3.is_var_yp0;
  EXPECT_TRUE(d3_yp0);
  auto d3_par = dae3.is_var_par;
  EXPECT_TRUE(d3_par);
  auto d3_fwd = dae3.use_fwd_sens;
  EXPECT_TRUE(d3_fwd);
}
}  // namespace idas_system_test
