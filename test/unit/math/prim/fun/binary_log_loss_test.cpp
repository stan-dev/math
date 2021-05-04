#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
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
}

TEST(MathFunctions, binary_log_loss_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::binary_log_loss(0, nan)));
  EXPECT_TRUE(std::isnan(stan::math::binary_log_loss(1, nan)));
}

TEST(MathFunctions, binary_log_loss_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::binary_log_loss;
    return binary_log_loss(x1, x2);
  };

  std::vector<int> std_in1{1, 3, 1};
  Eigen::VectorXd in2(3);
  in2 << 1.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, std_in1, in2);

  std::vector<std::vector<int>> std_std_in1{std_in1, std_in1, std_in1};
  Eigen::MatrixXd mat_in2 = in2.replicate(1, 3);
  stan::test::binary_scalar_tester(f, std_std_in1, mat_in2);
}
