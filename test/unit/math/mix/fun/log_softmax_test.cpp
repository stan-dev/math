#include <test/unit/math/test_ad.hpp>
#include <vector>

template <typename T>
void test_log_softmax(const T& x) {
  auto f = [](const auto& x) { return stan::math::log_softmax(x); };
  stan::test::expect_ad(f, x);
  stan::test::expect_ad_matvar(f, x);
}

TEST(MathMixMatFun, logSoftmax) {
  using Eigen::VectorXd;

  VectorXd x0(0);  // error case
  test_log_softmax(x0);

  VectorXd x1(1);
  x1 << 0;
  test_log_softmax(x1);

  VectorXd x2(2);
  x2 << -1, 1;
  test_log_softmax(x2);

  VectorXd x3(3);
  x3 << -1, 1, 10;
  test_log_softmax(x3);

  VectorXd x3b(3);
  x3b << 0, 1, 2;
  test_log_softmax(x3b);

  VectorXd x3c(3);
  x3c << 2, 1, 1;
  test_log_softmax(x3c);
}
