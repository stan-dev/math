#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsLkjCorr, fvar_double) {
  using stan::math::fvar;
  boost::random::mt19937 rng;
  int K = 4;
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, Eigen::Dynamic> Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  for (int i = 0; i < K * K; i++)
    Sigma(i).d_ = 1.0;
  fvar<double> eta = stan::math::uniform_rng(0, 2, rng);
  fvar<double> f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_, stan::math::lkj_corr_lpdf(Sigma, eta).val_);
  EXPECT_FLOAT_EQ(2.5177896, stan::math::lkj_corr_lpdf(Sigma, eta).d_);
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_, stan::math::lkj_corr_lpdf(Sigma, eta).val_);
  EXPECT_FLOAT_EQ(f.d_, stan::math::lkj_corr_lpdf(Sigma, eta).d_);
}

TEST(ProbDistributionsLkjCorrCholesky, fvar_double) {
  using stan::math::fvar;
  boost::random::mt19937 rng;
  int K = 4;
  Eigen::Matrix<fvar<double>, Eigen::Dynamic, Eigen::Dynamic> Sigma(K, K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  for (int i = 0; i < K * K; i++)
    Sigma(i).d_ = 1.0;
  fvar<double> eta = stan::math::uniform_rng(0, 2, rng);
  fvar<double> f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_, stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).val_);
  EXPECT_FLOAT_EQ(6.7766843, stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).d_);
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_, stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).val_);
  EXPECT_FLOAT_EQ(3, stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).d_);
}

TEST(ProbDistributionsLkjCorr, fvar_fvar_double) {
  using stan::math::fvar;
  boost::random::mt19937 rng;
  int K = 4;
  Eigen::Matrix<fvar<fvar<double> >, Eigen::Dynamic, Eigen::Dynamic> Sigma(K,
                                                                           K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  for (int i = 0; i < K * K; i++)
    Sigma(i).d_.val_ = 1.0;
  fvar<fvar<double> > eta = stan::math::uniform_rng(0, 2, rng);
  fvar<fvar<double> > f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_.val_, stan::math::lkj_corr_lpdf(Sigma, eta).val_.val_);
  EXPECT_FLOAT_EQ(2.5177896, stan::math::lkj_corr_lpdf(Sigma, eta).d_.val_);
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_.val_, stan::math::lkj_corr_lpdf(Sigma, eta).val_.val_);
  EXPECT_FLOAT_EQ(f.d_.val_, stan::math::lkj_corr_lpdf(Sigma, eta).d_.val_);
}

TEST(ProbDistributionsLkjCorrCholesky, fvar_fvar_double) {
  using stan::math::fvar;
  boost::random::mt19937 rng;
  int K = 4;
  Eigen::Matrix<fvar<fvar<double> >, Eigen::Dynamic, Eigen::Dynamic> Sigma(K,
                                                                           K);
  Sigma.setZero();
  Sigma.diagonal().setOnes();
  for (int i = 0; i < K * K; i++)
    Sigma(i).d_.val_ = 1.0;
  fvar<fvar<double> > eta = stan::math::uniform_rng(0, 2, rng);
  fvar<fvar<double> > f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_.val_,
                  stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).val_.val_);
  EXPECT_FLOAT_EQ(6.7766843,
                  stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).d_.val_);
  eta = 1.0;
  f = stan::math::do_lkj_constant(eta, K);
  EXPECT_FLOAT_EQ(f.val_.val_,
                  stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).val_.val_);
  EXPECT_FLOAT_EQ(3, stan::math::lkj_corr_cholesky_lpdf(Sigma, eta).d_.val_);
}
