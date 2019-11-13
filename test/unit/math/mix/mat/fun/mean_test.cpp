#include <test/unit/math/test_ad.hpp>

template <typename T>
void expect_mean(const T& m) {
  auto f = [](const auto& x) { return stan::math::mean(x); };
  Eigen::VectorXd v(m.size());
  Eigen::RowVectorXd rv(m.size());
  for (int i = 0; i < m.size(); ++i) {
    v(i) = m(i);
    rv(i) = m(i);
  }
  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
  stan::test::expect_ad(f, m);
}

TEST(MathMixMatFun, mean) {
  Eigen::MatrixXd a(0, 0);
  expect_mean(a);

  Eigen::MatrixXd b(1, 1);
  b << 12;
  expect_mean(b);

  Eigen::MatrixXd c(3, 1);
  c << 100, 0, -3;
  expect_mean(c);

  Eigen::MatrixXd d(1, 3);
  d << 100, 0, -3;
  expect_mean(d);

  Eigen::MatrixXd e(3, 2);
  e << -100, 0, 1, 20, -40, 2;
  expect_mean(e);
}
