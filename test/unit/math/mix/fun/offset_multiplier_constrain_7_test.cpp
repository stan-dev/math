#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, offset_multiplier_consistent_sizes) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    stan::return_type_t<decltype(x1), decltype(x2), decltype(x3)> lp = 0;
    return stan::math::offset_multiplier_constrain<false>(x1, x2, x3, lp);
  };

  double xd = 1.0;
  Eigen::VectorXd x(4);
  x << 1.0, 2.0, 3.0, 4.0;

  double mud = 2.0;
  Eigen::RowVectorXd mu(4);
  mu << 1.0, 2.0, 3.0, 4.0;

  double sigmad = 1.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << 1.0, 2.0, 3.0, 4.0;

  stan::test::expect_ad(f, xd, mu, sigma);
  stan::test::expect_ad(f, x, mud, sigma);
  stan::test::expect_ad(f, x, mu, sigmad);

  stan::test::expect_ad_matvar(f, xd, mu, sigma);
  stan::test::expect_ad_matvar(f, x, mud, sigma);
  stan::test::expect_ad_matvar(f, x, mu, sigmad);
}
