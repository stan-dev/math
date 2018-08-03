#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixTrace, fd) {
  using stan::math::fvar;
  using stan::math::matrix_fd;
  using stan::math::trace;

  matrix_fd a(2, 2);
  a << -1.0, 2.0, 5.0, 10.0;
  a(0, 0).d_ = 1.0;
  a(0, 1).d_ = 1.0;
  a(1, 0).d_ = 1.0;
  a(1, 1).d_ = 1.0;

  fvar<double> s = trace(a);
  EXPECT_FLOAT_EQ(9.0, s.val_);
  EXPECT_FLOAT_EQ(2.0, s.d_);
}
TEST(AgradFwdMatrixTrace, ffd) {
  using stan::math::fvar;
  using stan::math::matrix_ffd;
  using stan::math::trace;

  matrix_ffd a(2, 2);
  a << -1.0, 2.0, 5.0, 10.0;
  a(0, 0).d_ = 1.0;
  a(0, 1).d_ = 1.0;
  a(1, 0).d_ = 1.0;
  a(1, 1).d_ = 1.0;

  fvar<fvar<double> > s = trace(a);
  EXPECT_FLOAT_EQ(9.0, s.val_.val());
  EXPECT_FLOAT_EQ(2.0, s.d_.val());
}
