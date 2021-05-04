#include <test/unit/math/test_ad.hpp>
#include <vector>

template <typename T>
void expect_tail(const T& x, int n) {
  auto f = [](int i) {
    return [=](const auto& y) { return stan::math::tail(y, i); };
  };
  Eigen::VectorXd v = stan::test::to_vector(x);
  Eigen::RowVectorXd rv = stan::test::to_row_vector(x);
  stan::test::expect_ad(f(n), x);
  stan::test::expect_ad(f(n), v);
  stan::test::expect_ad(f(n), rv);
}

TEST(MathMixMatFun, tail) {
  std::vector<double> a{};
  expect_tail(a, 0);
  expect_tail(a, 1);

  std::vector<double> b{1};
  expect_tail(b, 0);
  expect_tail(b, 1);
  expect_tail(b, 2);

  std::vector<double> v{1, 2, 3};
  for (int n = 0; n < 5; ++n) {
    expect_tail(v, n);
  }
}

TEST(MathMixMatFun, tailEig) {
  Eigen::VectorXd a(0);
  expect_tail(a, 0);
  expect_tail(a, 1);

  Eigen::VectorXd b(1);
  b << 1;
  expect_tail(b, 0);
  expect_tail(b, 1);
  expect_tail(b, 2);

  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  for (int n = 0; n < 5; ++n) {
    expect_tail(v, n);
  }
}
