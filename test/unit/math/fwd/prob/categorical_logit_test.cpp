#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsCategoricalLogit, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::log_softmax;
  Matrix<fvar<double>, Dynamic, 1> theta(3, 1);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;
  Matrix<fvar<double>, Dynamic, 1> theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_,
                  stan::math::categorical_logit_lpmf(1, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_,
                  stan::math::categorical_logit_lpmf(2, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_,
                  stan::math::categorical_logit_lpmf(3, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_,
                  stan::math::categorical_logit_lpmf(1, theta).d_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_,
                  stan::math::categorical_logit_lpmf(2, theta).d_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_,
                  stan::math::categorical_logit_lpmf(3, theta).d_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_double_vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::log_softmax;
  Matrix<fvar<double>, Dynamic, 1> theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_lpmf(ns, theta).val_);

  Matrix<fvar<double>, Dynamic, 1> theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_ + theta_log_softmax[1].val_
                      + theta_log_softmax[0].val_,
                  stan::math::categorical_logit_lpmf(ms, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_ + theta_log_softmax[1].d_
                      + theta_log_softmax[0].d_,
                  stan::math::categorical_logit_lpmf(ms, theta).d_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::log_softmax;
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta(3, 1);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta_log_softmax
      = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_,
                  stan::math::categorical_logit_lpmf(1, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val_,
                  stan::math::categorical_logit_lpmf(2, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val_,
                  stan::math::categorical_logit_lpmf(3, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_,
                  stan::math::categorical_logit_lpmf(1, theta).d_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val_,
                  stan::math::categorical_logit_lpmf(2, theta).d_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val_,
                  stan::math::categorical_logit_lpmf(3, theta).d_.val_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_double_vectorized) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::log_softmax;
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_lpmf(ns, theta).val_.val_);

  Matrix<fvar<fvar<double> >, Dynamic, 1> theta_log_softmax
      = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_
                      + theta_log_softmax[1].val_.val_
                      + theta_log_softmax[0].val_.val_,
                  stan::math::categorical_logit_lpmf(ms, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_ + theta_log_softmax[1].d_.val_
                      + theta_log_softmax[0].d_.val_,
                  stan::math::categorical_logit_lpmf(ms, theta).d_.val_);
}
