#include <test/unit/math/test_ad.hpp>
#include <tuple>

TEST(MathMixFun, interp1PchipEvalIdentityPlaceholder) {
  Eigen::VectorXd xk(3);
  Eigen::VectorXd yk(3);
  xk << 0.0, 1.0, 3.0;
  yk << 2.0, 4.0, 8.0;
  auto W = stan::math::interp1_pchip_setup(xk, yk);

  auto f = [W](const auto& x) { return stan::math::interp1_pchip_eval(x, W); };
  stan::test::expect_ad(f, -1.5);
  stan::test::expect_ad(f, 0.0);
  stan::test::expect_ad(f, 1.25);
  stan::test::expect_ad(f, 7.0);
}
