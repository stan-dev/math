#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(test_unit_math_test_ad, to_std_vector) {
  Eigen::VectorXd u(0);
  EXPECT_EQ(0, stan::test::to_std_vector(u).size());

  Eigen::VectorXd x(3);
  x << 1, 2, 3;
  std::vector<double> y = stan::test::to_std_vector(x);
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(x(i), y[i]);
}

TEST(test_unit_math_test_ad, to_eigen_vector) {
  std::vector<double> u;
  EXPECT_EQ(0, stan::test::to_eigen_vector(u).size());

  std::vector<double> v{1, 2, 3};
  Eigen::VectorXd vv = stan::test::to_eigen_vector(v);
  EXPECT_EQ(3, vv.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(v[i], vv(i));
}

TEST(test_unit_math_test_ad, test_ad_unary) {
  Eigen::MatrixXd x(2, 2);
  x << 1.9, 0.3, 0.3, 1.7;

  auto g = [](const auto& u) { return stan::math::inverse(u); };
  stan::test::expect_ad(g, x);

  // need functor because can't infer template param F of expect_ad
  // stan::test::expect_ad(inverse, x);  // won't compile
}

TEST(test_unit_math_test_ad, test_ad_binary) {
  auto g
      = [](const auto& u, const auto& v) { return stan::math::multiply(u, v); };

  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  double y = -3;
  stan::test::expect_ad(g, x, y);

  double y2 = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(g, x, y2);
}
