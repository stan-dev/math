#include <test/unit/math/test_ad.hpp>
#include <tuple>

TEST(MathMixFun, interp1CubicHermiteEvalIdentityPlaceholder) {
  Eigen::VectorXd xk(3);
  Eigen::VectorXd yk(3);
  Eigen::VectorXd dydxk(3);
  xk << 0.0, 1.0, 3.0;
  yk << 2.0, 4.0, 8.0;
  dydxk << 1.0, 2.0, 1.5;
  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);

  auto f = [W](const auto& x) { return stan::math::interp1_cubic_hermite_eval(x, W); };
  stan::test::expect_ad(f, -1.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 1.25);
  stan::test::expect_ad(f, 7.0);
}
