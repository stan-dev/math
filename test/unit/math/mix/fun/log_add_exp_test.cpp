#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>
#include <iostream>

TEST(MathMixMatFun, logAddExp) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::log_add_exp(x, y);
  };
  // Test with finite values
  Eigen::VectorXd x1(2);
  x1 << 2.0, 1.0;
  Eigen::VectorXd y1(2);
  y1 << 3.0, 2.0;
  stan::test::expect_ad(f, x1, y1);

  // Test with negative infinity
  
  stan::test::expect_ad(f, stan::math::NEGATIVE_INFTY, 1.0);
  stan::test::expect_ad(f, 1.0, stan::math::NEGATIVE_INFTY);

  // Test with infinity
  stan::test::expect_ad(f, stan::math::INFTY, stan::math::INFTY);
}

TEST(MathMixMatFun, log_add_exp_elementwise_values) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::log_add_exp(x, y);
  };

  Eigen::VectorXd x1(2);
  x1 << 2.0, 1.0;
  Eigen::VectorXd y1(2);
  y1 << 3.0, 2.0;
  stan::test::expect_ad(f, x1, y1);

  Eigen::VectorXd x2(2);
  x2 << 0.5, -1.0;
  Eigen::VectorXd y2(2);
  y2 << 1.0, 2.0;
  stan::test::expect_ad(f, x2, y2);

  // Test with infinity
  Eigen::VectorXd x3(2);
  x3 << std::numeric_limits<double>::infinity(), 1.0;
  Eigen::VectorXd y3(2);
  y3 << 2.0, std::numeric_limits<double>::infinity();
  stan::test::expect_ad(f, x3, y3);
}

TEST(MathMixMatFun, log_add_exp_edge_cases) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::log_add_exp(x, y);
  };

  stan::test::expect_ad(f, stan::math::NEGATIVE_INFTY, 1.0);
  stan::test::expect_ad(f, 1.0, stan::math::NEGATIVE_INFTY);
  stan::test::expect_ad(f, stan::math::INFTY, stan::math::INFTY);
}

TEST(MathMixMatFun, log_add_exp_mismatched_sizes) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::log_add_exp(x, y);
  };

  std::vector<double> x{1.0, 2.0};
  std::vector<double> y{1.0, 2.0, 3.0};
  
  stan::test::expect_ad(f, x, y);
  stan::test::expect_ad(f, y, x);
}
