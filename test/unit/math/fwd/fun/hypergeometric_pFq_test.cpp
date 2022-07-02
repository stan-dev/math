#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(primScalFun, grad_2F2_fd) {
  using stan::math::fvar;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;
  using stan::math::vector_fd;

  vector_fd fd_a(2);
  vector_d d_a(2);
  fd_a.val() << 4, 2;
  d_a << 4, 2;
  fd_a.d() << 1, 1;

  vector_fd fd_b(2);
  vector_d d_b(2);
  fd_b << 6, 3;
  d_b << 6, 3;
  fd_b.d() << 1, 1;

  fvar<double> fd_z = fvar<double>(4, 1);
  double d_z = 4;

  double a_adj = 3.924636646666071 + 6.897245961898751;
  double b_adj = -2.775051002566842 - 4.980095849781222;
  double z_adj = 4.916522138006060;

  // fvar, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_a, fd_b, fd_z).d_,
                  a_adj + b_adj + z_adj);

  // fvar, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_a, fd_b, d_z).d_, a_adj + b_adj);

  // fvar, double, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_a, d_b, d_z).d_, a_adj);

  // fvar, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_a, d_b, fd_z).d_, a_adj + z_adj);

  // double, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, fd_b, fd_z).d_, b_adj + z_adj);

  // double, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, fd_b, d_z).d_, b_adj);

  // double, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, d_b, fd_z).d_, z_adj);
}

TEST(primScalFun, grad_2F2_ffd) {
  using stan::math::fvar;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;
  using stan::math::vector_ffd;

  vector_ffd ffd_a(2);
  vector_d d_a(2);
  ffd_a.val() << 4, 2;
  d_a << 4, 2;
  ffd_a.val().d() << 1, 1;

  vector_ffd ffd_b(2);
  vector_d d_b(2);
  ffd_b.val() << 6, 3;
  d_b << 6, 3;
  ffd_b.val().d() << 1, 1;

  fvar<fvar<double>> ffd_z;
  ffd_z.val_.val_ = 4;
  ffd_z.val_.d_ = 1;
  double d_z = 4;

  double a_adj = 3.924636646666071 + 6.897245961898751;
  double b_adj = -2.775051002566842 - 4.980095849781222;
  double z_adj = 4.916522138006060;

  // fvar, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffd_a, ffd_b, ffd_z).val_.d_,
                  a_adj + b_adj + z_adj);

  // fvar, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffd_a, ffd_b, d_z).val_.d_, a_adj + b_adj);

  // fvar, double, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffd_a, d_b, d_z).val_.d_, a_adj);

  // fvar, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffd_a, d_b, ffd_z).val_.d_, a_adj + z_adj);

  // double, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, ffd_b, ffd_z).val_.d_, b_adj + z_adj);

  // double, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, ffd_b, d_z).val_.d_, b_adj);

  // double, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, d_b, ffd_z).val_.d_, z_adj);
}
