#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(MathMixMatFun, divide) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::divide(x, y); };

  double x1 = 10;
  double x2 = -2;
  stan::test::expect_ad(f, x1, x2);

  Eigen::VectorXd v(1);
  v << 10;
  stan::test::expect_ad(f, v, x2);

  Eigen::RowVectorXd rv = v;
  stan::test::expect_ad(f, rv, x2);

  Eigen::MatrixXd m(1, 1);
  m << 10;
  stan::test::expect_ad(f, m, x2);

  Eigen::MatrixXd p(3, 2);
  p << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, p, x2);

  Eigen::VectorXd w(3);
  w << 100, 0, -3;
  stan::test::expect_ad(f, w, x2);

  Eigen::VectorXd u(4);
  u << 100, 0, -3, 4;
  stan::test::expect_ad(f, u, x2);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, v0, x1);
  stan::test::expect_ad(f, rv0, x1);
  stan::test::expect_ad(f, m00, x1);

  stan::test::expect_ad(f, u, 0.0);
  stan::test::expect_ad(f, rv, 0.0);
  stan::test::expect_ad(f, p, 0.0);

  Eigen::RowVectorXd rv4(4);
  rv4 << -5, 10, 7, 8.2;
  stan::test::expect_ad(f, rv4, x2);

  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (double value : {inf, -inf, nan}) {
    stan::test::expect_ad(f, u, value);
    stan::test::expect_ad(f, rv4, value);
    stan::test::expect_ad(f, p, value);
    stan::math::vector_d v1 = u;
    stan::math::row_vector_d rv1 = rv4;
    stan::math::matrix_d m1 = p;
    v1(1) = value;
    rv1(1) = value;
    m1(1, 1) = value;
    // XXX TODO: These fail, because there is no check_finite usage, and the
    // nans/infs propagate through autodiff and finite diff differently.
    // In particular:
    // f[x,y] = x/y
    // f_x[x,y] = 1/y
    // f_x[nan,1] = 1 (autodiff answer)
    // f_x[nan,1] ~ (1/n) (f[nan+1/n,1] - f[nan,1])
    //            = (1/n) (nan - nan)
    //            = nan (finite diff answer)
    // stan::test::expect_ad(f, v1, 2.0);
    // stan::test::expect_ad(f, rv1, 2.0);
    // stan::test::expect_ad(f, m1, 2.0);
  }
}
