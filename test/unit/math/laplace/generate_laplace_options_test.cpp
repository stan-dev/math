#include <stan/math/mix/functor/laplace_marginal_density_estimator.hpp>
#include <gtest/gtest.h>

namespace stan::math::test {

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromSizeDefaults) {
  constexpr int theta_0_size = 4;
  const auto defaults = laplace_options_default{};

  const auto options = stan::math::generate_laplace_options(theta_0_size);
  const auto& theta_0 = std::get<0>(options);

  EXPECT_EQ(theta_0_size, theta_0.size());
  EXPECT_TRUE(theta_0.isApprox(Eigen::VectorXd::Zero(theta_0_size)));
  EXPECT_DOUBLE_EQ(defaults.tolerance, std::get<1>(options));
  EXPECT_EQ(defaults.max_num_steps, std::get<2>(options));
  EXPECT_EQ(defaults.solver, std::get<3>(options));
  EXPECT_EQ(defaults.line_search.max_iterations, std::get<4>(options));
  EXPECT_EQ(static_cast<int>(defaults.allow_fallthrough), std::get<5>(options));
}

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromSizeZero) {
  const auto options = stan::math::generate_laplace_options(0);
  const auto& theta_0 = std::get<0>(options);

  EXPECT_EQ(0, theta_0.size());
}

TEST(LaplaceMarginalDensityEstimator, GenerateLaplaceOptionsFromThetaDefaults) {
  Eigen::VectorXd theta_0(3);
  theta_0 << -1.2, 0.5, 2.3;
  const auto defaults = laplace_options_default{};

  const auto options = stan::math::generate_laplace_options(theta_0);

  EXPECT_TRUE(theta_0.isApprox(std::get<0>(options)));
  EXPECT_DOUBLE_EQ(defaults.tolerance, std::get<1>(options));
  EXPECT_EQ(defaults.max_num_steps, std::get<2>(options));
  EXPECT_EQ(defaults.solver, std::get<3>(options));
  EXPECT_EQ(defaults.line_search.max_iterations, std::get<4>(options));
  EXPECT_EQ(static_cast<int>(defaults.allow_fallthrough), std::get<5>(options));
}

}  // namespace stan::math::test
