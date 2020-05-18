#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, binary_log_loss) {
  EXPECT_FLOAT_EQ(0.0, stan::math::binary_log_loss(0, 0.0));
  EXPECT_FLOAT_EQ(0.0, stan::math::binary_log_loss(1, 1.0));
  EXPECT_FLOAT_EQ(-log(0.5), stan::math::binary_log_loss(0, 0.5));
  EXPECT_FLOAT_EQ(-log(0.5), stan::math::binary_log_loss(1, 0.5));
  EXPECT_FLOAT_EQ(-log(0.75), stan::math::binary_log_loss(0, 0.25));
  EXPECT_FLOAT_EQ(-log(0.75), stan::math::binary_log_loss(1, 0.75));

  Eigen::VectorXi in1(3);
  in1 << 0, 0, 0;
  Eigen::VectorXd in2(3);
  in2 << 0.5, 0.5, 0.5;
  Eigen::MatrixXd in_mat(3,3);
  in_mat.fill(0.5);

  Eigen::VectorXd out(3);
  out = stan::math::binary_log_loss(in1, 0.5);
  out = stan::math::binary_log_loss(0, in2);
  out = stan::math::binary_log_loss(in1, in2);
  out = stan::math::binary_log_loss(in1, in_mat.diagonal());

  std::vector<int> st_in1{0, 0, 0};
  std::vector<std::vector<int>> stst_in1{st_in1, st_in1, st_in1};
  std::vector<double> st_in2{0.5, 0.5, 0.5};
  std::vector<std::vector<double>> stst_in2{st_in2, st_in2, st_in2};
  std::vector<double> st_out(3);
  std::vector<std::vector<double>> stst_out(3);
  std::vector<Eigen::VectorXi> st_eig1{in1, in1, in1};
  std::vector<Eigen::VectorXd> st_eig2{in2, in2, in2};
  std::vector<Eigen::VectorXd> st_eigout(3);
  st_out = stan::math::binary_log_loss(st_in1, 0.5);
  stst_out = stan::math::binary_log_loss(stst_in1, 0.5);
  st_out = stan::math::binary_log_loss(0, st_in2);
  stst_out = stan::math::binary_log_loss(0, stst_in2);
  stst_out = stan::math::binary_log_loss(stst_in1, stst_in2);
  st_eigout = stan::math::binary_log_loss(st_eig1, 0.5);
  st_eigout = stan::math::binary_log_loss(st_eig1, st_eig2);
  st_eigout = stan::math::binary_log_loss(0, st_eig2);
  double res = -log(0.5);
  for(int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(res, out[i]);
    EXPECT_FLOAT_EQ(res, st_out[i]);
    EXPECT_FLOAT_EQ(res, stst_out[i][i]);
    EXPECT_FLOAT_EQ(res, st_eigout[i][i]);
  }
}

TEST(MathFunctions, binary_log_loss_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::binary_log_loss(0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::binary_log_loss(1, nan)));
}
