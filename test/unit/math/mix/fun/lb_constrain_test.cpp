#include <test/unit/math/test_ad.hpp>
#include <limits>

void expect_lb_constrain(double x, double lb) {
  auto f1 = [](const auto& x, const auto& lb) {
    return stan::math::lb_constrain(x, lb);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain(x, lb, lp);
    return lp;
  };

  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
}

TEST(mathMixScalFun, lbConstrain) {
  expect_lb_constrain(-1, 2);
  expect_lb_constrain(2, 4);
}

TEST(mathMixMatFun, lb_mat_constrain) {
  using stan::scalar_type_t;
  using stan::math::lb_constrain;
  using stan::math::promote_scalar_t;

  auto f1 = [](const auto& x, const auto& lb) {
    return stan::math::lb_constrain(x, lb);
  };
  auto f2 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    return stan::math::lb_constrain(x, lb, lp);
  };
  auto f3 = [](const auto& x, const auto& lb) {
    stan::return_type_t<decltype(x), decltype(lb)> lp = 0;
    stan::math::lb_constrain(x, lb, lp);
    return lp;
  };

  Eigen::MatrixXd A(2, 2);
  A << -5.0, 2.0, 4.0, -2.0;
  Eigen::MatrixXd lbm(2, 2);
  lbm << 1.0, -5.0, 0.0, 100.0;

  double lbd = -5.0;

  stan::test::expect_ad(f1, A, lbm);
  stan::test::expect_ad(f1, A, lbd);
  stan::test::expect_ad_matvar(f1, A, lbm);
  stan::test::expect_ad_matvar(f1, A, lbd);
  stan::test::expect_ad(f2, A, lbm);
  stan::test::expect_ad(f2, A, lbd);
  stan::test::expect_ad_matvar(f2, A, lbm);
  stan::test::expect_ad_matvar(f2, A, lbd);
  stan::test::expect_ad(f3, A, lbm);
  stan::test::expect_ad(f3, A, lbd);
  stan::test::expect_ad_matvar(f3, A, lbm);
  stan::test::expect_ad_matvar(f3, A, lbd);
}
