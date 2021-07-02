#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(MathMixMatFun, divide_vec_scal) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::divide(x, y); };

  double x1 = 10;
  double x2 = -2;

  Eigen::VectorXd v(1);
  v << 10;
  stan::test::expect_ad(f, x2, v);
  stan::test::expect_ad_matvar(f, x2, v);

  Eigen::RowVectorXd rv = v;
  stan::test::expect_ad(f, x2, rv);
  stan::test::expect_ad_matvar(f, x2, rv);

  Eigen::MatrixXd m(1, 1);
  m << 10;
  stan::test::expect_ad(f, x2, m);
  stan::test::expect_ad_matvar(f, x2, m);

  Eigen::MatrixXd p(3, 2);
  p << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, x2, p);
  stan::test::expect_ad_matvar(f, x2, m);
  Eigen::VectorXd w(3);
  w << 100, 2, -3;
  stan::test::expect_ad(f, x2, w);
  stan::test::expect_ad_matvar(f, x2, m);

  Eigen::VectorXd u(4);
  u << 100, 0.5, -3, 4;
  stan::test::expect_ad(f, x2, u);
  stan::test::expect_ad_matvar(f, x2, m);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, x1, v0);
  stan::test::expect_ad(f, x1, rv0);
  stan::test::expect_ad(f, x1, m00);

  stan::test::expect_ad_matvar(f, x1, v0);
  stan::test::expect_ad_matvar(f, x1, rv0);
  stan::test::expect_ad_matvar(f, x1, m00);

  Eigen::RowVectorXd rv4(4);
  rv4 << -5, 10, 7, 8.2;
  stan::test::expect_ad(f, x2, rv4);
  stan::test::expect_ad_matvar(f, x2, rv4);

  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (double value : {inf, -inf, nan}) {
    stan::test::expect_ad(f, value, u);
    stan::test::expect_ad(f, value, rv4);
    stan::test::expect_ad(f, value, p);

    stan::test::expect_ad_matvar(f, value, u);
    stan::test::expect_ad_matvar(f, value, rv4);
    stan::test::expect_ad_matvar(f, value, p);
  }
}
