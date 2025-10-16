#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <vector>

TEST(StanMathRev, vd) {
  using stan::math::var;
  harm_osc_ode_data_fun harm_osc;
  stan::math::internal::integrate_ode_std_vector_interface_adapter<
      harm_osc_ode_data_fun>
      harm_osc_adapted(harm_osc);

  std::vector<var> theta = {0.15};
  std::vector<double> y = {1.0, 0.5};

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  double t = 1.0;

  Eigen::Matrix<var, Eigen::Dynamic, 1> out1
      = stan::math::to_vector(harm_osc(t, y, theta, x, x_int, nullptr));
  Eigen::Matrix<var, Eigen::Dynamic, 1> out2
      = harm_osc_adapted(t, stan::math::to_vector(y), nullptr, theta, x, x_int);

  stan::math::sum(out1).grad();
  Eigen::VectorXd adjs1(theta.size());
  for (size_t i = 0; i < theta.size(); ++i)
    adjs1(i) = theta[i].adj();

  stan::math::set_zero_all_adjoints();

  stan::math::sum(out2).grad();
  Eigen::VectorXd adjs2(theta.size());
  for (size_t i = 0; i < theta.size(); ++i)
    adjs2(i) = theta[i].adj();

  EXPECT_MATRIX_FLOAT_EQ(out1.val(), out2.val());
  EXPECT_MATRIX_FLOAT_EQ(adjs1, adjs2);
}

TEST(StanMathRev, dv) {
  using stan::math::var;
  harm_osc_ode_data_fun harm_osc;
  stan::math::internal::integrate_ode_std_vector_interface_adapter<
      harm_osc_ode_data_fun>
      harm_osc_adapted(harm_osc);

  std::vector<double> theta = {0.15};
  std::vector<var> y = {1.0, 0.5};

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  double t = 1.0;

  Eigen::Matrix<var, Eigen::Dynamic, 1> out1
      = stan::math::to_vector(harm_osc(t, y, theta, x, x_int, nullptr));
  Eigen::Matrix<var, Eigen::Dynamic, 1> out2
      = harm_osc_adapted(t, stan::math::to_vector(y), nullptr, theta, x, x_int);

  stan::math::sum(out1).grad();
  Eigen::VectorXd adjs1(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    adjs1(i) = y[i].adj();

  stan::math::set_zero_all_adjoints();

  stan::math::sum(out2).grad();
  Eigen::VectorXd adjs2(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    adjs2(i) = y[i].adj();

  EXPECT_MATRIX_FLOAT_EQ(out1.val(), out2.val());
  EXPECT_MATRIX_FLOAT_EQ(adjs1, adjs2);
}

TEST(StanMathRev, vv) {
  using stan::math::var;
  harm_osc_ode_data_fun harm_osc;
  stan::math::internal::integrate_ode_std_vector_interface_adapter<
      harm_osc_ode_data_fun>
      harm_osc_adapted(harm_osc);

  std::vector<var> theta = {0.15};
  std::vector<var> y = {1.0, 0.5};

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  double t = 1.0;

  Eigen::Matrix<var, Eigen::Dynamic, 1> out1
      = stan::math::to_vector(harm_osc(t, y, theta, x, x_int, nullptr));
  Eigen::Matrix<var, Eigen::Dynamic, 1> out2
      = harm_osc_adapted(t, stan::math::to_vector(y), nullptr, theta, x, x_int);

  stan::math::sum(out1).grad();
  Eigen::VectorXd adjs_theta_1(theta.size());
  for (size_t i = 0; i < theta.size(); ++i)
    adjs_theta_1(i) = theta[i].adj();
  Eigen::VectorXd adjs_y_1(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    adjs_y_1(i) = y[i].adj();

  stan::math::set_zero_all_adjoints();

  stan::math::sum(out2).grad();
  Eigen::VectorXd adjs_theta_2(theta.size());
  for (size_t i = 0; i < theta.size(); ++i)
    adjs_theta_2(i) = theta[i].adj();
  Eigen::VectorXd adjs_y_2(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    adjs_y_2(i) = y[i].adj();

  EXPECT_MATRIX_FLOAT_EQ(out1.val(), out2.val());
  EXPECT_MATRIX_FLOAT_EQ(adjs_theta_1, adjs_theta_2);
  EXPECT_MATRIX_FLOAT_EQ(adjs_y_1, adjs_y_2);
}
