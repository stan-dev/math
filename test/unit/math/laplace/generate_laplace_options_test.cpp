#include <stan/math/prim/fun/generate_laplace_options.hpp>
#include <gtest/gtest.h>

namespace stan::math::test {

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromSizeDefaults) {
  constexpr int theta_0_size = 4;
  const auto options = stan::math::generate_laplace_options(theta_0_size);
  const auto& theta_0 = std::get<0>(options);

  EXPECT_EQ(theta_0_size, theta_0.size());
  EXPECT_TRUE(theta_0.isApprox(Eigen::VectorXd::Zero(theta_0_size)));
  EXPECT_DOUBLE_EQ(1.49012e-08, std::get<1>(options));
  EXPECT_EQ(500, std::get<2>(options));
  EXPECT_EQ(1, std::get<3>(options));
  EXPECT_EQ(1000, std::get<4>(options));
  EXPECT_EQ(1, std::get<5>(options));
}

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromSizeZero) {
  const auto options = stan::math::generate_laplace_options(0);
  const auto& theta_0 = std::get<0>(options);

  EXPECT_EQ(0, theta_0.size());
}

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromThetaDefaults) {
  Eigen::VectorXd theta_0(3);
  theta_0 << -1.2, 0.5, 2.3;
  const auto options = stan::math::generate_laplace_options(theta_0);

  EXPECT_TRUE(theta_0.isApprox(std::get<0>(options)));
  EXPECT_DOUBLE_EQ(1.49012e-08, std::get<1>(options));
  EXPECT_EQ(500, std::get<2>(options));
  EXPECT_EQ(1, std::get<3>(options));
  EXPECT_EQ(1000, std::get<4>(options));
  EXPECT_EQ(1, std::get<5>(options));
}

}  // namespace stan::math::test
