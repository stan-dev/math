#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <limits>


using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::log_softmax;
using stan::math::log_softmax;

TEST(ProbDistributionsCategoricalLogit,fvar_var) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<var>,Dynamic,1> theta(3,1);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;
  Matrix<fvar<var>,Dynamic,1> theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val(), 
                  stan::math::categorical_logit_log(1,theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val(), 
                  stan::math::categorical_logit_log(2,theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val(), 
                  stan::math::categorical_logit_log(3,theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val(), 
                  stan::math::categorical_logit_log(1,theta).d_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val(), 
                  stan::math::categorical_logit_log(2,theta).d_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val(), 
                  stan::math::categorical_logit_log(3,theta).d_.val());
}

TEST(ProbDistributionsCategoricalLogit,fvar_var_vectorized) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<var>,Dynamic,1> theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns,theta).val_.val());

  Matrix<fvar<var>,Dynamic,1> theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val() + theta_log_softmax[1].val_.val() 
                  + theta_log_softmax[0].val_.val(),
                  stan::math::categorical_logit_log(ms,theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val() + theta_log_softmax[1].d_.val() 
                  + theta_log_softmax[0].d_.val(),
                  stan::math::categorical_logit_log(ms,theta).d_.val());
}

TEST(ProbDistributionsCategoricalLogit,fvar_fvar_var) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<fvar<var> >,Dynamic,1> theta(3,1);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;
  Matrix<fvar<fvar<var> >,Dynamic,1> theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_.val(), 
                  stan::math::categorical_logit_log(1,theta).val_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val_.val(), 
                  stan::math::categorical_logit_log(2,theta).val_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val_.val(), 
                  stan::math::categorical_logit_log(3,theta).val_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_.val(), 
                  stan::math::categorical_logit_log(1,theta).d_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val_.val(),
                  stan::math::categorical_logit_log(2,theta).d_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val_.val(),
                  stan::math::categorical_logit_log(3,theta).d_.val_.val());
}

TEST(ProbDistributionsCategoricalLogit,fvar_fvar_var_vectorized) {
  using stan::math::var;
  using stan::math::fvar;
  Matrix<fvar<fvar<var> >,Dynamic,1> theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns,theta).val_.val_.val());

  Matrix<fvar<fvar<var> >,Dynamic,1> theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_.val() 
                  + theta_log_softmax[1].val_.val_.val() 
                  + theta_log_softmax[0].val_.val_.val(),
                  stan::math::categorical_logit_log(ms,theta).val_.val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_.val() 
                  + theta_log_softmax[1].d_.val_.val() 
                  + theta_log_softmax[0].d_.val_.val(),
                  stan::math::categorical_logit_log(ms,theta).d_.val_.val());
}
