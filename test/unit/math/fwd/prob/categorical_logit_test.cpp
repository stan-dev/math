#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsCategoricalLogit, fvar_double) {
  using stan::math::vector_fd;
  using stan::math::log_softmax;
  vector_fd theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;
  vector_fd theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_,
                  stan::math::categorical_logit_log(1, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_,
                  stan::math::categorical_logit_log(2, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_,
                  stan::math::categorical_logit_log(3, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_,
                  stan::math::categorical_logit_log(1, theta).d_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_,
                  stan::math::categorical_logit_log(2, theta).d_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_,
                  stan::math::categorical_logit_log(3, theta).d_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_double_vectorized) {
  using stan::math::vector_fd;
  using stan::math::log_softmax;
  vector_fd theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns, theta).val_);

  vector_fd theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_ + theta_log_softmax[1].val_
                      + theta_log_softmax[0].val_,
                  stan::math::categorical_logit_log(ms, theta).val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_ + theta_log_softmax[1].d_
                      + theta_log_softmax[0].d_,
                  stan::math::categorical_logit_log(ms, theta).d_);

  std::vector<vector_fd> theta_arr(3);
  theta_arr[0] = (vector_fd(3) << 1, 4, 3).finished();
  theta_arr[1] = (vector_fd(3) << -1, -3, 2.1).finished();
  theta_arr[2] = (vector_fd(3) << 0.3, 0.4, -2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_ = i;
    theta_arr[1][i].d_ = i * 3;
    theta_arr[2][i].d_ = i / 2;
  }

  std::vector<vector_fd> theta_arr_log_softmax = log_softmax(theta_arr);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].val_
                    + theta_arr_log_softmax[1][1].val_
                    + theta_arr_log_softmax[2][0].val_,
                  stan::math::categorical_logit_log(ms, theta_arr).val_);
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].d_
                    + theta_arr_log_softmax[1][1].d_
                    + theta_arr_log_softmax[2][0].d_,
                  stan::math::categorical_logit_log(ms, theta_arr).d_);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].val_
                    + theta_arr_log_softmax[1][1].val_
                    + theta_arr_log_softmax[2][1].val_,
                  stan::math::categorical_logit_log(2, theta_arr).val_);
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].d_
                    + theta_arr_log_softmax[1][1].d_
                    + theta_arr_log_softmax[2][1].d_,
                  stan::math::categorical_logit_log(2, theta_arr).d_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_double) {
  using stan::math::vector_ffd;
  using stan::math::log_softmax;
  vector_ffd theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;
  vector_ffd theta_log_softmax
      = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_,
                  stan::math::categorical_logit_log(1, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val_,
                  stan::math::categorical_logit_log(2, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val_,
                  stan::math::categorical_logit_log(3, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_,
                  stan::math::categorical_logit_log(1, theta).d_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val_,
                  stan::math::categorical_logit_log(2, theta).d_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val_,
                  stan::math::categorical_logit_log(3, theta).d_.val_);
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_double_vectorized) {
  using stan::math::vector_ffd;
  using stan::math::log_softmax;
  vector_ffd theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns, theta).val_.val_);

  vector_ffd theta_log_softmax
      = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val_
                      + theta_log_softmax[1].val_.val_
                      + theta_log_softmax[0].val_.val_,
                  stan::math::categorical_logit_log(ms, theta).val_.val_);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val_ + theta_log_softmax[1].d_.val_
                      + theta_log_softmax[0].d_.val_,
                  stan::math::categorical_logit_log(ms, theta).d_.val_);

  std::vector<vector_ffd> theta_arr(3);
  theta_arr[0] = (vector_ffd(3) << 1, 4, 3).finished();
  theta_arr[1] = (vector_ffd(3) << -1, -3, 2.1).finished();
  theta_arr[2] = (vector_ffd(3) << 0.3, 0.4, -2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_.val_ = i;
    theta_arr[1][i].d_.val_ = i * 3;
    theta_arr[2][i].d_.val_ = i / 2;
  }

  std::vector<vector_ffd> theta_arr_log_softmax = log_softmax(theta_arr);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].val_.val_
                      + theta_arr_log_softmax[1][1].val_.val_
                      + theta_arr_log_softmax[2][0].val_.val_,
                  stan::math::categorical_logit_log(ms, theta_arr).val_.val_);
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].d_.val_
                      + theta_arr_log_softmax[1][1].d_.val_
                      + theta_arr_log_softmax[2][0].d_.val_,
                  stan::math::categorical_logit_log(ms, theta_arr).d_.val_);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].val_.val_
                      + theta_arr_log_softmax[1][1].val_.val_
                      + theta_arr_log_softmax[2][1].val_.val_,
                  stan::math::categorical_logit_log(2, theta_arr).val_.val_);
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].d_.val_
                      + theta_arr_log_softmax[1][1].d_.val_
                      + theta_arr_log_softmax[2][1].d_.val_,
                  stan::math::categorical_logit_log(2, theta_arr).d_.val_);
}
