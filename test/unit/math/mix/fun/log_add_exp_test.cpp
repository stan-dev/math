#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

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

  Eigen::VectorXd result = f(x3, y3);
  EXPECT_TRUE(std::isinf(result[0]));  // Expect infinity for the first element
  EXPECT_TRUE(std::isinf(result[1]));
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

TEST(MathMixMatFun, log_add_exp_container_tests) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::log_add_exp(x, y);
  };

  // Test with Eigen::MatrixXd
  Eigen::MatrixXd x_row(1, 2);
  x_row << 2.0, 1.0;
  Eigen::MatrixXd y_row(1, 2);
  y_row << 3.0, 2.0;

  stan::test::expect_ad(f, x_row, y_row);

  // Additional tests with mismatched sizes
  Eigen::MatrixXd x_mismatch(2, 1);
  x_mismatch << 0.5, -1.0;
  Eigen::MatrixXd y_mismatch(1, 3);
  y_mismatch << 1.0, 2.0, 3.0;

  EXPECT_THROW(stan::math::log_add_exp(x_mismatch, y_mismatch),
               std::invalid_argument);
}
