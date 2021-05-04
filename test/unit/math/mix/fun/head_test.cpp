#include <test/unit/math/test_ad.hpp>
#include <vector>

template <typename T>
void expect_head(const T& x, int n) {
  auto f = [](int i) {
    return [=](const auto& y) { return stan::math::head(y, i); };
  };
  Eigen::VectorXd v = stan::test::to_vector(x);
  Eigen::RowVectorXd rv = stan::test::to_row_vector(x);
  stan::test::expect_ad(f(n), x);
  stan::test::expect_ad(f(n), v);
  stan::test::expect_ad(f(n), rv);
}

TEST(MathMixMatFun, head) {
  std::vector<double> a{};
  expect_head(a, 0);
  expect_head(a, 1);

  std::vector<double> b{1};
  expect_head(b, 0);
  expect_head(b, 1);
  expect_head(b, 2);

  std::vector<double> v{1, 2, 3};
  for (int n = 0; n < 5; ++n) {
    expect_head(v, n);
  }
}
