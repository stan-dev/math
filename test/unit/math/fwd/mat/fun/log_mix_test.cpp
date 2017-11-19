#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::fvar;
using stan::math::log_mix;

TEST(AgradFwdMatrixLogMix, vector_fd) {
  using stan::math::vector_fd;

  vector_fd dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  vector_fd prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<double> a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_);
}

TEST(AgradFwdMatrixLogMix, row_vector_fd) {
  using stan::math::row_vector_fd;

  row_vector_fd dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  row_vector_fd prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<double> a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_);
}

TEST(AgradFwdMatrixLogMix, std_vector_fd) {
  std::vector<fvar<double> > dens(4);
  dens[0].val_ = -1.0;
  dens[1].val_ = -2.0;
  dens[2].val_ = -3.0;
  dens[3].val_ = -4.0;
  dens[0].d_ = 1.0;
  dens[1].d_ = 1.0;
  dens[2].d_ = 1.0;
  dens[3].d_ = 1.0;

  std::vector<fvar<double> > prob(4);
  prob[0].val_ = 0.15;
  prob[1].val_ = 0.70;
  prob[2].val_ = 0.10;
  prob[3].val_ = 0.05;
  prob[0].d_ = 1.0;
  prob[1].d_ = 1.0;
  prob[2].d_ = 1.0;
  prob[3].d_ = 1.0;

  fvar<double> a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_);
}

TEST(AgradFwdMatrixLogMix, vector_ffd) {
  using stan::math::vector_ffd;

  vector_ffd dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  vector_ffd prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<fvar<double> > a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_.val_);
}

TEST(AgradFwdMatrixLogMix, row_vector_ffd) {
  using stan::math::row_vector_ffd;

  row_vector_ffd dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  row_vector_ffd prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<fvar<double> > a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_.val_);
}

TEST(AgradFwdMatrixLogMix, std_vector_ffd) {
  std::vector<fvar<fvar<double> > > dens(4);
  dens[0].val_ = -1.0;
  dens[1].val_ = -2.0;
  dens[2].val_ = -3.0;
  dens[3].val_ = -4.0;
  dens[0].d_ = 1.0;
  dens[1].d_ = 1.0;
  dens[2].d_ = 1.0;
  dens[3].d_ = 1.0;

  std::vector<fvar<fvar<double> > > prob(4);
  prob[0].val_ = 0.15;
  prob[1].val_ = 0.70;
  prob[2].val_ = 0.10;
  prob[3].val_ = 0.05;
  prob[0].d_ = 1.0;
  prob[1].d_ = 1.0;
  prob[2].d_ = 1.0;
  prob[3].d_ = 1.0;

  fvar<fvar<double> > a = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(-1.85911088, a.val_.val_);
  EXPECT_FLOAT_EQ(4.66673118, a.d_.val_);
}
