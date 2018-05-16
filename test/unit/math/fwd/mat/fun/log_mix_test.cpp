#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::fvar;
using stan::math::log_mix;
using stan::math::row_vector_fd;
using stan::math::row_vector_ffd;
using stan::math::vector_fd;
using stan::math::vector_ffd;

template <typename T_a, typename T_b>
void fvar_test(T_a a, T_b b) {
  a[0].val_ = 0.235;
  a[1].val_ = 0.152;
  a[2].val_ = 0.359;
  a[3].val_ = 0.254;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b[0].val_ = -5.581;
  b[1].val_ = -6.254;
  b[2].val_ = -3.987;
  b[3].val_ = -10.221;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;

  fvar<double> out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_, -4.8474302);
  EXPECT_FLOAT_EQ(out.d_, 4.093988932);

  T_b b2(4), b3(4);

  b2[0].val_ = -6.785;
  b2[1].val_ = -4.351;
  b2[2].val_ = -5.847;
  b2[3].val_ = -7.362;
  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;

  b3[0].val_ = -7.251;
  b3[1].val_ = -10.510;
  b3[2].val_ = -12.302;
  b3[3].val_ = -3.587;
  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;

  std::vector<T_b> c{b, b2, b3};

  fvar<double> std_out = log_mix(a, c);

  EXPECT_FLOAT_EQ(std_out.val_, -15.457609);
  EXPECT_FLOAT_EQ(std_out.d_, 15.164879);
}

template <typename T_a, typename T_b>
void fvarfvar_test(T_a a, T_b b) {
  a[0].val_ = 0.235;
  a[1].val_ = 0.152;
  a[2].val_ = 0.359;
  a[3].val_ = 0.254;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b[0].val_ = -5.581;
  b[1].val_ = -6.254;
  b[2].val_ = -3.987;
  b[3].val_ = -10.221;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;

  fvar<fvar<double> > out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val_, -4.8474302);
  EXPECT_FLOAT_EQ(out.d_.val_, 4.093988932);

  T_b b2(4), b3(4);

  b2[0].val_ = -6.785;
  b2[1].val_ = -4.351;
  b2[2].val_ = -5.847;
  b2[3].val_ = -7.362;
  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;

  b3[0].val_ = -7.251;
  b3[1].val_ = -10.510;
  b3[2].val_ = -12.302;
  b3[3].val_ = -3.587;
  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;

  std::vector<T_b> c{b, b2, b3};

  fvar<fvar<double> > std_out = log_mix(a, c);

  EXPECT_FLOAT_EQ(std_out.val_.val_, -15.457609);
  EXPECT_FLOAT_EQ(std_out.d_.val_, 15.164879);
}

TEST(AgradFwdMatrixLogMix, fdValues) {
  /**
   * Test that all possible combinations of inputs return
   * the same result and derivatives.
   */

  vector_fd vecd_prob(4);
  vector_fd vecd_dens(4);
  row_vector_fd row_vecd_prob(4);
  row_vector_fd row_vecd_dens(4);
  std::vector<fvar<double> > std_prob(4);
  std::vector<fvar<double> > std_dens(4);

  fvar_test(vecd_prob, vecd_dens);
  fvar_test(vecd_prob, row_vecd_dens);
  fvar_test(vecd_prob, std_dens);

  fvar_test(row_vecd_prob, vecd_dens);
  fvar_test(row_vecd_prob, row_vecd_dens);
  fvar_test(row_vecd_prob, std_dens);

  fvar_test(std_prob, vecd_dens);
  fvar_test(std_prob, row_vecd_dens);
  fvar_test(std_prob, std_dens);
}

TEST(AgradFwdMatrixLogMix, ffdValues) {
  /**
   * Test that all possible combinations of inputs return
   * the same result and derivatives.
   */

  vector_ffd vecd_prob(4);
  vector_ffd vecd_dens(4);
  row_vector_ffd row_vecd_prob(4);
  row_vector_ffd row_vecd_dens(4);
  std::vector<fvar<fvar<double> > > std_prob(4);
  std::vector<fvar<fvar<double> > > std_dens(4);

  fvarfvar_test(vecd_prob, vecd_dens);
  fvarfvar_test(vecd_prob, row_vecd_dens);
  fvarfvar_test(vecd_prob, std_dens);

  fvarfvar_test(row_vecd_prob, vecd_dens);
  fvarfvar_test(row_vecd_prob, row_vecd_dens);
  fvarfvar_test(row_vecd_prob, std_dens);

  fvarfvar_test(std_prob, vecd_dens);
  fvarfvar_test(std_prob, row_vecd_dens);
  fvarfvar_test(std_prob, std_dens);
}

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
