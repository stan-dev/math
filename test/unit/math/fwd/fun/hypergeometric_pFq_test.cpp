#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(primScalFun, grad_2F2_fd) {
  using stan::math::vector_fd;
  using stan::math::vector_d;
  using stan::math::fvar;
  using stan::math::hypergeometric_pFq;

  vector_fd fd_p(2);
  vector_d d_p(2);
  fd_p.val() << 4, 2;
  d_p << 4, 2;
  fd_p.d() << 1, 1;

  vector_fd fd_q(2);
  vector_d d_q(2);
  fd_q << 6, 3;
  d_q << 6, 3;
  fd_q.d() << 1, 1;

  fvar<double> fd_z = fvar<double>(4, 1);
  double d_z = 4;

  double p_adj = 3.924636646666071 + 6.897245961898751;
  double q_adj = -2.775051002566842 - 4.980095849781222;
  double z_adj = 4.916522138006060;

  // fvar, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_p, fd_q, fd_z).d_,
                  p_adj + q_adj + z_adj);

  // fvar, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_p, fd_q, d_z).d_,
                  p_adj + q_adj);

  // fvar, double, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_p, d_q, d_z).d_,
                  p_adj);

  // fvar, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(fd_p, d_q, fd_z).d_,
                  p_adj + z_adj);

  // double, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, fd_q, fd_z).d_,
                  q_adj + z_adj);

  // double, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, fd_q, d_z).d_, q_adj);

  // double, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, d_q, fd_z).d_, z_adj);
}

TEST(primScalFun, grad_2F2_ffd) {
  using stan::math::vector_ffd;
  using stan::math::vector_fd;
  using stan::math::fvar;
  using stan::math::hypergeometric_pFq;

  vector_ffd fd_p(2);
  vector_fd d_p(2);
  fd_p.val() << 4, 2;
  d_p.val() << 4, 2;
  fd_p.val().d() << 1, 1;

  vector_ffd fd_q(2);
  vector_fd d_q(2);
  fd_q.val() << 6, 3;
  d_q.val() << 6, 3;
  fd_q.val().d() << 1, 1;

  fvar<fvar<double>> ffd_z;;
  ffd_z.val_.val_ = 4;
  ffd_z.val_.d_ = 1;
  double d_z = 4;

  // fvar, fvar, fvar
  fvar<fvar<double>> result = hypergeometric_pFq(fd_p, fd_q, ffd_z);

  double adj = 3.924636646666071 + 6.897245961898751
                 -2.775051002566842 - 4.980095849781222
                 + 4.916522138006060;

  EXPECT_FLOAT_EQ(result.val_.d_, adj);
}
