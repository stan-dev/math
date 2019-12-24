#include <test/unit/math/test_ad.hpp>
#include <vector>

void expect_cumulative_sum(std::vector<double>& x) {
  auto f = [](const auto& y) { return stan::math::cumulative_sum(y); };
  Eigen::VectorXd v = Eigen::Map<Eigen::VectorXd>(x.data(), x.size());
  Eigen::RowVectorXd rv = Eigen::Map<Eigen::RowVectorXd>(x.data(), x.size());
  stan::test::expect_ad(f, x);
  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
}

TEST(MathMixMatFun, cumulativeSum) {
  std::vector<double> a;
  expect_cumulative_sum(a);

  std::vector<double> b{1.7};
  expect_cumulative_sum(b);

  std::vector<double> c{5.9, -1.2};
  expect_cumulative_sum(c);

  std::vector<double> d{5.9, -1.2, 192.13456};
  expect_cumulative_sum(d);
}
