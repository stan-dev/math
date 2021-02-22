#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, lub_constrain_scalar) {
  auto f1 = [](const auto& x, const auto& lb, const auto& ub) {
    return stan::math::lub_constrain(x, lb, ub);
  };
  auto f2 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    return stan::math::lub_constrain(x, lb, ub, lp);
  };
  auto f3 = [](const auto& x, const auto& lb, const auto& ub) {
    stan::return_type_t<decltype(x), decltype(lb), decltype(ub)> lp = 0;
    stan::math::lub_constrain(x, lb, ub, lp);
    return lp;
  };

  double x1 = 0.7;
  double x2 = -1.1;
  double lb = -2.0;
  double ub = 1.9;

  stan::test::expect_ad(f1, x1, lb, ub);
  stan::test::expect_ad(f1, x2, lb, ub);

  stan::test::expect_ad(f2, x1, lb, ub);
  stan::test::expect_ad(f2, x2, lb, ub);

  stan::test::expect_ad(f3, x1, lb, ub);
  stan::test::expect_ad(f3, x2, lb, ub);

  // Error cases
  stan::test::expect_ad(f1, x1, lb, lb);
  stan::test::expect_ad(f1, x2, lb, lb);

  stan::test::expect_ad(f2, x1, lb, lb);
  stan::test::expect_ad(f2, x2, lb, lb);

  stan::test::expect_ad(f3, x1, lb, lb);
  stan::test::expect_ad(f3, x2, lb, lb);
}
