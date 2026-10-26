#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, scaleMultiply) {
  auto f = [](const auto& s, const auto& a, const auto& b) {
    return stan::math::scale_multiply(s, a, b);
  };
  double s = 0.7;
  Eigen::MatrixXd a(3, 2);
  a << 1, -2, 0.5, 3, -1.5, 2;
  Eigen::MatrixXd b(2, 2);
  b << 0.25, -1, 2, 0.5;
  Eigen::VectorXd v(2);
  v << -0.5, 1.5;
  Eigen::RowVectorXd rv(3);
  rv << 1, 0.5, -2;
  // reverse mode only: scale_multiply has no fvar overload yet
  constexpr bool reverse_only = true;
  stan::test::expect_ad<reverse_only>(f, s, a, b);
  stan::test::expect_ad<reverse_only>(f, s, a, v);
  stan::test::expect_ad<reverse_only>(f, s, rv, a);
  Eigen::MatrixXd empty(0, 2);
  stan::test::expect_ad<reverse_only>(f, s, empty, b);
}
