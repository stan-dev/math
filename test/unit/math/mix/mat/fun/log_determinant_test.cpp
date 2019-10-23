#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, logDeterminant) {
  auto f = [](const auto& x) { return stan::math::log_determinant(x); };

  // repeat symmetric pos def tests
  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);

  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (const auto& y : stan::test::ar_test_cov_matrices(1, 3, rho)) {
      stan::test::expect_ad(f, y);
    }
  }

  // asymmetric tests
  Eigen::MatrixXd d(2, 3);
  d << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, d);

  Eigen::MatrixXd e(2, 2);
  e << 0, 1, 2, 3;
  stan::test::expect_ad(f, e);

  Eigen::MatrixXd g(2, 2);
  e << 2, 3, 5, 7;
  stan::test::expect_ad(f, g);
}

#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(AgradMixMatrixLogDeterminant, fv_1stDeriv) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(0.0, 1.0);
  fvar<var> b(1.0, 2.0);
  fvar<var> c(2.0, 2.0);
  fvar<var> d(3.0, 2.0);

  matrix_fv v(2, 2);
  v << a, b, c, d;

  fvar<var> det;
  det = log_determinant(v);
  EXPECT_FLOAT_EQ(std::log(2.0), det.val_.val());
  EXPECT_FLOAT_EQ(1.5, det.d_.val());

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  det.val_.grad(q, h);
  EXPECT_FLOAT_EQ(-1.5, h[0]);
  EXPECT_FLOAT_EQ(1, h[1]);
  EXPECT_FLOAT_EQ(.5, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, fv_2ndDeriv) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_fv;
  using stan::math::var;

  fvar<var> a(0.0, 1.0);
  fvar<var> b(1.0, 2.0);
  fvar<var> c(2.0, 2.0);
  fvar<var> d(3.0, 2.0);
  matrix_fv v(2, 2);
  v << a, b, c, d;

  fvar<var> det;
  det = log_determinant(v);

  AVEC q = createAVEC(a.val(), b.val(), c.val(), d.val());
  VEC h;
  det.d_.grad(q, h);
  EXPECT_FLOAT_EQ(1.25, h[0]);
  EXPECT_FLOAT_EQ(-.5, h[1]);
  EXPECT_FLOAT_EQ(0.25, h[2]);
  EXPECT_FLOAT_EQ(-.5, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, fv_exception) {
  using stan::math::log_determinant;
  using stan::math::matrix_fv;

  EXPECT_THROW(log_determinant(matrix_fv(2, 3)), std::invalid_argument);
}
TEST(AgradMixMatrixLogDeterminant, ffv_1stDeriv) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(0.0, 1.0);
  fvar<fvar<var> > b(1.0, 2.0);
  fvar<fvar<var> > c(2.0, 2.0);
  fvar<fvar<var> > d(3.0, 2.0);

  matrix_ffv v(2, 2);
  v << a, b, c, d;

  fvar<fvar<var> > det;
  det = log_determinant(v);
  EXPECT_FLOAT_EQ(std::log(2.0), det.val_.val().val());
  EXPECT_FLOAT_EQ(1.5, det.d_.val().val());

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  det.val_.val().grad(q, h);
  EXPECT_FLOAT_EQ(-1.5, h[0]);
  EXPECT_FLOAT_EQ(1, h[1]);
  EXPECT_FLOAT_EQ(.5, h[2]);
  EXPECT_FLOAT_EQ(0.0, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, ffv_2ndDeriv_1) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(0.0, 1.0);
  fvar<fvar<var> > b(1.0, 2.0);
  fvar<fvar<var> > c(2.0, 2.0);
  fvar<fvar<var> > d(3.0, 2.0);
  matrix_ffv v(2, 2);
  v << a, b, c, d;

  fvar<fvar<var> > det;
  det = log_determinant(v);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  det.val().d_.grad(q, h);
  EXPECT_FLOAT_EQ(0, h[0]);
  EXPECT_FLOAT_EQ(0, h[1]);
  EXPECT_FLOAT_EQ(0, h[2]);
  EXPECT_FLOAT_EQ(0, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, ffv_2ndDeriv_2) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(0.0, 1.0);
  fvar<fvar<var> > b(1.0, 2.0);
  fvar<fvar<var> > c(2.0, 2.0);
  fvar<fvar<var> > d(3.0, 2.0);
  matrix_ffv v(2, 2);
  v << a, b, c, d;

  fvar<fvar<var> > det;
  det = log_determinant(v);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  det.d_.val().grad(q, h);
  EXPECT_FLOAT_EQ(1.25, h[0]);
  EXPECT_FLOAT_EQ(-.5, h[1]);
  EXPECT_FLOAT_EQ(0.25, h[2]);
  EXPECT_FLOAT_EQ(-.5, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, ffv_3rdDeriv) {
  using stan::math::fvar;
  using stan::math::log_determinant;
  using stan::math::matrix_ffv;
  using stan::math::var;

  fvar<fvar<var> > a(0.0, 1.0);
  fvar<fvar<var> > b(1.0, 2.0);
  fvar<fvar<var> > c(2.0, 2.0);
  fvar<fvar<var> > d(3.0, 2.0);
  a.val_.d_ = 1.0;
  b.val_.d_ = 1.0;
  c.val_.d_ = 1.0;
  d.val_.d_ = 1.0;

  matrix_ffv v(2, 2);
  v << a, b, c, d;

  fvar<fvar<var> > det;
  det = log_determinant(v);

  AVEC q
      = createAVEC(a.val().val(), b.val().val(), c.val().val(), d.val().val());
  VEC h;
  det.d_.d_.grad(q, h);
  EXPECT_FLOAT_EQ(1.5, h[0]);
  EXPECT_FLOAT_EQ(-1.25, h[1]);
  EXPECT_FLOAT_EQ(-1, h[2]);
  EXPECT_FLOAT_EQ(0.75, h[3]);
}
TEST(AgradMixMatrixLogDeterminant, ffv_exception) {
  using stan::math::log_determinant;
  using stan::math::matrix_ffv;

  EXPECT_THROW(log_determinant(matrix_ffv(2, 3)), std::invalid_argument);
}
