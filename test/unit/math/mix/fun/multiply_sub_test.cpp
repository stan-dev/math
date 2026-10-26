#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, multiplySub) {
  auto f = [](const auto& a, const auto& b, const auto& c) {
    return stan::math::multiply_sub(a, b, c);
  };
  Eigen::MatrixXd a(3, 2);
  a << 1, -2, 0.5, 3, -1.5, 2;
  Eigen::MatrixXd b(2, 2);
  b << 0.25, -1, 2, 0.5;
  Eigen::MatrixXd c(3, 2);
  c << 0.1, 0.2, -0.3, 0.4, 0.5, -0.6;
  Eigen::VectorXd v(2);
  v << -0.5, 1.5;
  Eigen::VectorXd cv(3);
  cv << 1, 2, 3;
  // reverse mode only: multiply_sub has no fvar overload yet
  constexpr bool reverse_only = true;
  stan::test::expect_ad<reverse_only>(f, a, b, c);
  stan::test::expect_ad<reverse_only>(f, a, v, cv);
}
