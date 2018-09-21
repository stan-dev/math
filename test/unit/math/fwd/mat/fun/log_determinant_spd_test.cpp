#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixLogDeterminantSPD, fd) {
  using stan::math::fvar;
  using stan::math::log_determinant_spd;
  using stan::math::matrix_fd;

  matrix_fd v(2, 2);
  v << 3, 0, 0, 4;
  v(0, 0).d_ = 1.0;
  v(0, 1).d_ = 2.0;
  v(1, 0).d_ = 2.0;
  v(1, 1).d_ = 2.0;

  fvar<double> det;
  det = log_determinant_spd(v);
  EXPECT_FLOAT_EQ(std::log(12.0), det.val_);
  EXPECT_FLOAT_EQ(0.83333333, det.d_);
}

TEST(AgradFwdMatrixLogDeterminantSPD, fd_exception) {
  using stan::math::log_determinant_spd;
  using stan::math::matrix_fd;

  EXPECT_THROW(log_determinant_spd(matrix_fd(2, 3)), std::invalid_argument);
}
TEST(AgradFwdMatrixLogDeterminantSPD, ffd) {
  using stan::math::fvar;
  using stan::math::log_determinant_spd;
  using stan::math::matrix_ffd;

  fvar<fvar<double> > a, b, c, d;
  a.val_.val_ = 3.0;
  a.d_.val_ = 1.0;
  b.val_.val_ = 0.0;
  b.d_.val_ = 2.0;
  c.val_.val_ = 0.0;
  c.d_.val_ = 2.0;
  d.val_.val_ = 4.0;
  d.d_.val_ = 2.0;

  matrix_ffd v(2, 2);
  v << a, b, c, d;

  fvar<fvar<double> > det;
  det = log_determinant_spd(v);
  EXPECT_FLOAT_EQ(std::log(12.0), det.val_.val());
  EXPECT_FLOAT_EQ(0.83333333, det.d_.val());
}

TEST(AgradFwdMatrixLogDeterminantSPD, ffd_exception) {
  using stan::math::log_determinant_spd;
  using stan::math::matrix_ffd;

  EXPECT_THROW(log_determinant_spd(matrix_ffd(2, 3)), std::invalid_argument);
}
