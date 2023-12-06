#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, inverseLDLTtest) {
  auto f = [](const auto& x) {
    auto x_sym = stan::math::multiply(0.5, x + x.transpose());
    auto ldlt = stan::math::make_ldlt_factor(x_sym);
    return stan::math::inverse_ldlt(ldlt);
  };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd a(1, 1);
  a << 2;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd c(2, 2);
  c << 2, 3, 3, 7;
  stan::test::expect_ad(f, c);

  Eigen::MatrixXd h(3, 3);  // not positive definite
  h << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  stan::test::expect_ad(f, h);
}