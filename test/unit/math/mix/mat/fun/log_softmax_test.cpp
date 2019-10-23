#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, logSoftmax) {
  auto f = [](const auto& x) { return stan::math::log_softmax(x); };
  Eigen::VectorXd x0(0);  // error case
  stan::test::expect_ad(f, x0);

  Eigen::VectorXd x1(1);
  x1 << 0;
  stan::test::expect_ad(f, x1);

  Eigen::VectorXd x2(2);
  x2 << -1, 1;
  stan::test::expect_ad(f, x2);

  Eigen::VectorXd x3(3);
  x3 << -1, 1, 10;
  stan::test::expect_ad(f, x3);

  Eigen::VectorXd x3b(3);
  x3 << 0, 1, 2;
  stan::test::expect_ad(f, x3b);

  Eigen::VectorXd x3c(3);
  x3 << 2, 1, 1;
  stan::test::expect_ad(f, x3c);
}
