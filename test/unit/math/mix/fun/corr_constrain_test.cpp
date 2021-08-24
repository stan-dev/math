#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, corr_constrain) {
  auto f = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    return stan::math::corr_constrain<false>(x1, lp);
  };

  std::vector<double> x0 = {-1.0, 2.0, 3.0};
  Eigen::VectorXd x1(3);
  x1 << -1.0, -2.0, 0.3;
  Eigen::RowVectorXd x2(3);
  x2 << -1.0, -2.0, 0.3;
  stan::test::expect_ad(f, x0);
  stan::test::expect_ad(f, x1);
  stan::test::expect_ad(f, x2);
  stan::test::expect_ad_matvar(f, x1);
  stan::test::expect_ad_matvar(f, x2);

  std::vector<double> x4;
  Eigen::VectorXd x5(0);
  Eigen::RowVectorXd x6(0);

  stan::test::expect_ad(f, x4);
  stan::test::expect_ad(f, x5);
  stan::test::expect_ad(f, x6);

  stan::test::expect_ad_matvar(f, x5);
  stan::test::expect_ad_matvar(f, x6);
}

TEST(mathMixMatFun, corr_constrain_lp) {
  auto f1 = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    return stan::math::corr_constrain<true>(x1, lp);
  };

  auto f2 = [](const auto& x1) {
    stan::scalar_type_t<decltype(x1)> lp = 0.0;
    stan::math::corr_constrain<true>(x1, lp);
    return lp;
  };

  std::vector<double> x0 = {-1.0, 2.0, 3.0};
  Eigen::VectorXd x1(3);
  x1 << -1.0, -2.0, 0.3;
  Eigen::RowVectorXd x2(3);
  x2 << -1.0, -2.0, 0.3;
  stan::test::expect_ad(f1, x0);
  stan::test::expect_ad(f1, x1);
  stan::test::expect_ad(f1, x2);
  stan::test::expect_ad_matvar(f1, x1);
  stan::test::expect_ad_matvar(f1, x2);

  stan::test::expect_ad(f2, x0);
  stan::test::expect_ad(f2, x1);
  stan::test::expect_ad(f2, x2);
  stan::test::expect_ad_matvar(f2, x1);
  stan::test::expect_ad_matvar(f2, x2);

  std::vector<double> x4;
  Eigen::VectorXd x5(0);
  Eigen::RowVectorXd x6(0);

  stan::test::expect_ad(f1, x4);
  stan::test::expect_ad(f1, x5);
  stan::test::expect_ad(f1, x6);
  stan::test::expect_ad_matvar(f1, x5);
  stan::test::expect_ad_matvar(f1, x6);

  stan::test::expect_ad(f2, x4);
  stan::test::expect_ad(f2, x5);
  stan::test::expect_ad(f2, x6);
  stan::test::expect_ad_matvar(f2, x5);
  stan::test::expect_ad_matvar(f2, x6);
}
