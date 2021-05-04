#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, reverse_vector) {
  auto f = [](const auto& x) { return stan::math::reverse(x); };

  // 0 x 0
  Eigen::VectorXd x0(0);
  stan::test::expect_ad(f, x0);

  // 1 x 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(f, x1);

  // 4 x 4
  Eigen::VectorXd x(4);
  x << 1, 3, 4, -5;
  stan::test::expect_ad(f, x);
}

TEST(MathMixMatFun, reverse_row_vector) {
  auto f = [](const auto& x) { return stan::math::reverse(x); };

  // 0 x 0
  Eigen::RowVectorXd x0(0);
  stan::test::expect_ad(f, x0);

  // 1 x 1
  Eigen::RowVectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(f, x1);

  // 4 x 4
  Eigen::RowVectorXd x(4);
  x << 1, 3, 4, -5;
  stan::test::expect_ad(f, x);
}

TEST(MathMixMatFun, reverse_array) {
  auto f = [](const auto& x) { return stan::math::reverse(x); };

  // 0 x 0
  std::vector<double> x0(0);
  stan::test::expect_ad(f, x0);

  // 1 x 1
  std::vector<double> x1{1};
  stan::test::expect_ad(f, x1);

  // 4 x 4
  std::vector<double> x{1, 3, 4, -5};
  stan::test::expect_ad(f, x);
}
