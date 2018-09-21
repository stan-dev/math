#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixQrR, fd) {
  using stan::math::matrix_d;
  using stan::math::matrix_fd;
  matrix_fd m0(0, 0);
  matrix_d m2(3, 2);
  matrix_fd m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;
  m2 << 1, 2, 3, 4, 5, 6;
  m1(0, 0).d_ = 1.0;
  m1(0, 1).d_ = 1.0;
  m1(1, 0).d_ = 1.0;
  m1(1, 1).d_ = 1.0;
  m1(2, 0).d_ = 1.0;
  m1(2, 1).d_ = 1.0;

  using stan::math::qr_thin_R;
  using stan::math::transpose;
  EXPECT_THROW(qr_thin_R(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_thin_R(m1));

  matrix_fd res = qr_thin_R(m1);
  matrix_d res2 = qr_thin_R(m2);

  for (int i = 0; i < res2.rows(); i++)
    for (int j = 0; j < res2.cols(); j++)
      EXPECT_FLOAT_EQ(res2(i, j), res(i, j).val_);

  EXPECT_FLOAT_EQ(1.5212777, res(0, 0).d_);
  EXPECT_FLOAT_EQ(1.6371845, res(0, 1).d_);
  EXPECT_FLOAT_EQ(0, res(1, 0).d_);
  EXPECT_FLOAT_EQ(-0.21293451, res(1, 1).d_);
  // EXPECT_FLOAT_EQ(0, res(2, 0).d_);
  // EXPECT_FLOAT_EQ(0, res(2, 1).d_);
}

TEST(AgradFwdMatrixQrR, ffd) {
  using stan::math::matrix_d;
  using stan::math::matrix_ffd;
  matrix_ffd m0(0, 0);
  matrix_d m2(3, 2);
  matrix_ffd m1(3, 2);
  m1 << 1, 2, 3, 4, 5, 6;
  m2 << 1, 2, 3, 4, 5, 6;
  m1(0, 0).d_ = 1.0;
  m1(0, 1).d_ = 1.0;
  m1(1, 0).d_ = 1.0;
  m1(1, 1).d_ = 1.0;
  m1(2, 0).d_ = 1.0;
  m1(2, 1).d_ = 1.0;

  using stan::math::qr_thin_R;
  using stan::math::transpose;
  EXPECT_THROW(qr_thin_R(m0), std::invalid_argument);
  EXPECT_NO_THROW(qr_thin_R(m1));

  matrix_ffd res = qr_thin_R(m1);
  matrix_d res2 = qr_thin_R(m2);

  for (int i = 0; i < res2.rows(); i++)
    for (int j = 0; j < res2.cols(); j++)
      EXPECT_FLOAT_EQ(res2(i, j), res(i, j).val_.val_);

  EXPECT_FLOAT_EQ(1.5212777, res(0, 0).d_.val_);
  EXPECT_FLOAT_EQ(1.6371845, res(0, 1).d_.val_);
  EXPECT_FLOAT_EQ(0, res(1, 0).d_.val_);
  EXPECT_FLOAT_EQ(-0.21293451, res(1, 1).d_.val_);
  // EXPECT_FLOAT_EQ(0, res(2, 0).d_.val_);
  // EXPECT_FLOAT_EQ(0, res(2, 1).d_.val_);
}
