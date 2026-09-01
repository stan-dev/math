#include <test/unit/math/test_ad.hpp>
#include <algorithm>

template <typename T>
inline void expect_max(const T& m) {
  auto f = [](const auto& x) { return stan::math::max(x); };
  Eigen::VectorXd v(m.size());
  Eigen::RowVectorXd rv(m.size());
  for (int i = 0; i < m.size(); ++i) {
    v(i) = m(i);
    rv(i) = m(i);
  }
  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
  stan::test::expect_ad(f, m);
  stan::test::expect_ad_matvar(f, v);
  stan::test::expect_ad_matvar(f, rv);
  stan::test::expect_ad_matvar(f, m);
}

TEST(MathMixMatFun, max) {
  Eigen::MatrixXd a(0, 0);
  expect_max(a);
  Eigen::MatrixXd b(1, 1);
  b << 12;
  expect_max(b);
  Eigen::MatrixXd c(3, 1);
  c << 100, 0, -3;
  expect_max(c);
  Eigen::MatrixXd d(1, 3);
  d << 100, 0, -3;
  expect_max(d);
  Eigen::MatrixXd e(3, 2);
  e << -100, 0, 1, 20, -40, 2;
  expect_max(e);
  Eigen::MatrixXd ties(2, 2);
  ties << 10.5, 1.0, 10.5, 10.5;
  expect_max(ties);
  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd o(3);
  o << -inf, 5.0, inf;
  expect_max(o);
  double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd n1(3);
  n1 << 1.0, nan, 2.0;
  expect_max(n1);
  Eigen::MatrixXd n2(2, 2);
  n2 << nan, nan, nan, nan;
  expect_max(n2);
}
TEST(MathMixMatFun, max_binary) {
  auto f = [](const auto& x, const auto& y) { return stan::math::max(x, y); };
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 2.0, 1.0);
  stan::test::expect_ad(f, 3.0, 3.0);
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, nan);
}
