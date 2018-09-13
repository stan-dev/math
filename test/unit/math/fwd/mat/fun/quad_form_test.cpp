#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

using stan::math::check_not_nan;
using stan::math::fvar;

TEST(AgradFwdMatrixQuadForm, quad_form_mat_fd) {
  using stan::math::matrix_fd;
  using stan::math::quad_form;

  matrix_fd ad(4, 4);
  matrix_fd bd(4, 2);

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  ad(0, 0).d_ = 1.0;
  ad(0, 1).d_ = 1.0;
  ad(0, 2).d_ = 1.0;
  ad(0, 3).d_ = 1.0;
  ad(1, 0).d_ = 1.0;
  ad(1, 1).d_ = 1.0;
  ad(1, 2).d_ = 1.0;
  ad(1, 3).d_ = 1.0;
  ad(2, 0).d_ = 1.0;
  ad(2, 1).d_ = 1.0;
  ad(2, 2).d_ = 1.0;
  ad(2, 3).d_ = 1.0;
  ad(3, 0).d_ = 1.0;
  ad(3, 1).d_ = 1.0;
  ad(3, 2).d_ = 1.0;
  ad(3, 3).d_ = 1.0;
  bd(0, 0).d_ = 1.0;
  bd(0, 1).d_ = 1.0;
  bd(1, 0).d_ = 1.0;
  bd(1, 1).d_ = 1.0;
  bd(2, 0).d_ = 1.0;
  bd(2, 1).d_ = 1.0;
  bd(3, 0).d_ = 1.0;
  bd(3, 1).d_ = 1.0;

  // fvar<double> - fvar<double>
  matrix_fd resd = quad_form(ad, bd);
  EXPECT_FLOAT_EQ(26033, resd(0, 0).val_);
  EXPECT_FLOAT_EQ(3456, resd(0, 1).val_);
  EXPECT_FLOAT_EQ(3396, resd(1, 0).val_);
  EXPECT_FLOAT_EQ(725, resd(1, 1).val_);
  EXPECT_FLOAT_EQ(15226, resd(0, 0).d_);
  EXPECT_FLOAT_EQ(3429, resd(0, 1).d_);
  EXPECT_FLOAT_EQ(4233, resd(1, 0).d_);
  EXPECT_FLOAT_EQ(900, resd(1, 1).d_);
}

TEST(AgradFwdMatrixQuadForm, quad_form_vec_fd) {
  using stan::math::matrix_fd;
  using stan::math::quad_form;
  using stan::math::vector_fd;

  matrix_fd ad(4, 4);
  vector_fd bd(4);
  fvar<double> res;

  bd << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  ad(0, 0).d_ = 1.0;
  ad(0, 1).d_ = 1.0;
  ad(0, 2).d_ = 1.0;
  ad(0, 3).d_ = 1.0;
  ad(1, 0).d_ = 1.0;
  ad(1, 1).d_ = 1.0;
  ad(1, 2).d_ = 1.0;
  ad(1, 3).d_ = 1.0;
  ad(2, 0).d_ = 1.0;
  ad(2, 1).d_ = 1.0;
  ad(2, 2).d_ = 1.0;
  ad(2, 3).d_ = 1.0;
  ad(3, 0).d_ = 1.0;
  ad(3, 1).d_ = 1.0;
  ad(3, 2).d_ = 1.0;
  ad(3, 3).d_ = 1.0;
  bd(0).d_ = 1.0;
  bd(1).d_ = 1.0;
  bd(2).d_ = 1.0;
  bd(3).d_ = 1.0;

  // fvar<double> - fvar<double>
  res = quad_form(ad, bd);
  EXPECT_FLOAT_EQ(26033, res.val_);
  EXPECT_FLOAT_EQ(15226, res.d_);
}

TEST(AgradFwdMatrixQuadForm, quad_form_mat_ffd) {
  using stan::math::matrix_ffd;
  using stan::math::quad_form;

  matrix_ffd ad(4, 4);
  matrix_ffd bd(4, 2);

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  ad(0, 0).d_ = 1.0;
  ad(0, 1).d_ = 1.0;
  ad(0, 2).d_ = 1.0;
  ad(0, 3).d_ = 1.0;
  ad(1, 0).d_ = 1.0;
  ad(1, 1).d_ = 1.0;
  ad(1, 2).d_ = 1.0;
  ad(1, 3).d_ = 1.0;
  ad(2, 0).d_ = 1.0;
  ad(2, 1).d_ = 1.0;
  ad(2, 2).d_ = 1.0;
  ad(2, 3).d_ = 1.0;
  ad(3, 0).d_ = 1.0;
  ad(3, 1).d_ = 1.0;
  ad(3, 2).d_ = 1.0;
  ad(3, 3).d_ = 1.0;
  bd(0, 0).d_ = 1.0;
  bd(0, 1).d_ = 1.0;
  bd(1, 0).d_ = 1.0;
  bd(1, 1).d_ = 1.0;
  bd(2, 0).d_ = 1.0;
  bd(2, 1).d_ = 1.0;
  bd(3, 0).d_ = 1.0;
  bd(3, 1).d_ = 1.0;

  // fvar<fvar<double> > - fvar<fvar<double> >
  matrix_ffd resd = quad_form(ad, bd);
  EXPECT_FLOAT_EQ(26033, resd(0, 0).val_.val_);
  EXPECT_FLOAT_EQ(3456, resd(0, 1).val_.val_);
  EXPECT_FLOAT_EQ(3396, resd(1, 0).val_.val_);
  EXPECT_FLOAT_EQ(725, resd(1, 1).val_.val_);
  EXPECT_FLOAT_EQ(15226, resd(0, 0).d_.val_);
  EXPECT_FLOAT_EQ(3429, resd(0, 1).d_.val_);
  EXPECT_FLOAT_EQ(4233, resd(1, 0).d_.val_);
  EXPECT_FLOAT_EQ(900, resd(1, 1).d_.val_);
}

TEST(AgradFwdMatrixQuadForm, quad_form_vec_ffd) {
  using stan::math::matrix_ffd;
  using stan::math::quad_form;
  using stan::math::vector_ffd;

  matrix_ffd ad(4, 4);
  vector_ffd bd(4);
  fvar<fvar<double> > res;

  bd << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  ad(0, 0).d_ = 1.0;
  ad(0, 1).d_ = 1.0;
  ad(0, 2).d_ = 1.0;
  ad(0, 3).d_ = 1.0;
  ad(1, 0).d_ = 1.0;
  ad(1, 1).d_ = 1.0;
  ad(1, 2).d_ = 1.0;
  ad(1, 3).d_ = 1.0;
  ad(2, 0).d_ = 1.0;
  ad(2, 1).d_ = 1.0;
  ad(2, 2).d_ = 1.0;
  ad(2, 3).d_ = 1.0;
  ad(3, 0).d_ = 1.0;
  ad(3, 1).d_ = 1.0;
  ad(3, 2).d_ = 1.0;
  ad(3, 3).d_ = 1.0;
  bd(0).d_ = 1.0;
  bd(1).d_ = 1.0;
  bd(2).d_ = 1.0;
  bd(3).d_ = 1.0;

  // fvar<fvar<double> > - fvar<fvar<double> >
  res = quad_form(ad, bd);
  EXPECT_FLOAT_EQ(26033, res.val_.val_);
  EXPECT_FLOAT_EQ(15226, res.d_.val_);
}
