#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, repMatrix) {
  // y is scalar
  auto f = [](int m, int n) {
    return [=](const auto& y) {
      return stan::math::rep_matrix<decltype(y)>(y, m, n);
    };
  };

  // y is row vector or column vector
  auto g = [](int k) {
    return [=](const auto& y) {
      return stan::math::rep_matrix<decltype(y)>(y, k);
    };
  };

  double y = 3;
  stan::test::expect_ad(f(0, 0), y);
  stan::test::expect_ad(f(1, 1), y);
  stan::test::expect_ad(f(2, 3), y);

  // illegal arguments---test throw
  stan::test::expect_ad(f(-2, -1), y);

  Eigen::VectorXd a(3);
  a << 3, 3, 3;
  stan::test::expect_ad(g(2), a);
  stan::test::expect_ad_matvar(g(2), a);

  Eigen::RowVectorXd b(2);
  b << 2, 2;
  stan::test::expect_ad(g(3), b);
  stan::test::expect_ad_matvar(g(3), b);
}
