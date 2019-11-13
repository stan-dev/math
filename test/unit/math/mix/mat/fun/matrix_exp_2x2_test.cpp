#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, matrixExp2x2) {
  auto f = [](const auto& x) { return stan::math::matrix_exp_2x2(x); };

  // example from Moler & Van Loan, 2003
  Eigen::MatrixXd x(2, 2);
  double a = -1;
  double b = -17;
  x << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;
  stan::test::expect_ad(f, x);

  Eigen::MatrixXd y(2, 2);
  y << 3, 0, 0, 4;
  stan::test::expect_ad(f, y);
}
