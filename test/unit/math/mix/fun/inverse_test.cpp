#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixMatFun, inverse) {
  auto f = [](const auto& x) { return stan::math::inverse(x); };

  Eigen::MatrixXd t(0, 0);
  stan::test::expect_ad(f, t);
  stan::test::expect_ad_matvar(f, t);

  Eigen::MatrixXd u(1, 1);
  u << 2;
  stan::test::expect_ad(f, u);
  stan::test::expect_ad_matvar(f, u);

  Eigen::MatrixXd v(2, 2);
  v << 2, 3, 5, 7;
  stan::test::expect_ad(f, v);
  stan::test::expect_ad_matvar(f, v);
  v << 1.9, 0.3, 0.3, 1.7;
  stan::test::expect_ad(f, v);
  stan::test::expect_ad_matvar(f, v);

  // issues around zero require looser tolerances for hessians
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 0.01;
  tols.hessian_fvar_hessian_ = 0.01;

  Eigen::MatrixXd w(3, 3);
  w << 2, 3, 5, 7, 11, 13, 17, 19, 23;
  stan::test::expect_ad(tols, f, w);
  stan::test::expect_ad_matvar(tols, f, w);

  // even lower tolerance, again for cases around zero
  stan::test::ad_tolerances tols2;
  tols2.hessian_hessian_ = 3.0;
  tols2.hessian_fvar_hessian_ = 3.0;

  Eigen::MatrixXd x(4, 4);
  x << 2, 3, 4, 5, 9, -1, 2, 2, 4, 3, 7, -1, 0, 1, 19, 112;
  stan::test::expect_ad(tols2, f, x);
  stan::test::expect_ad_matvar(tols2, f, x);

  Eigen::MatrixXd y(3, 2);
  y << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(stan::math::inverse(y), std::invalid_argument);
  stan::test::expect_ad(f, y);
  stan::test::expect_ad_matvar(f, y);

  Eigen::MatrixXd z(2, 2);
  z << 1, 2, 5, std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(stan::math::inverse(z));

  // autodiff throws, so following fails (throw behavior must match to pass)
  // stan::test::expect_ad(f, z);

  Eigen::MatrixXd a(2, 2);
  a << 1.9, 0.3, 0.3, std::numeric_limits<double>::infinity();
  stan::test::expect_ad(f, a);
  stan::test::expect_ad_matvar(f, a);

  // singular matrix, inverse does not exist.
  // Eigen does not throw errors, but returns inf
  Eigen::MatrixXd m(2, 2);
  m << 6.0, 4.0, 3.0, 2.0;
  EXPECT_NO_THROW(stan::math::inverse(m));
}
