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
