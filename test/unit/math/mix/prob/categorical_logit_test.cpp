#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsCategoricalLogit, fvar_var) {
  using stan::math::vector_fv;
  using stan::math::log_softmax;

  vector_fv theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;
  vector_fv theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val(),
                  stan::math::categorical_logit_log(1, theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val(),
                  stan::math::categorical_logit_log(2, theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val(),
                  stan::math::categorical_logit_log(3, theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val(),
                  stan::math::categorical_logit_log(1, theta).d_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val(),
                  stan::math::categorical_logit_log(2, theta).d_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val(),
                  stan::math::categorical_logit_log(3, theta).d_.val());
}

TEST(ProbDistributionsCategoricalLogit, fvar_var_vectorized) {
  using stan::math::vector_fv;
  using stan::math::log_softmax;

  vector_fv theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0, stan::math::categorical_logit_log(ns, theta).val_.val());

  vector_fv theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val()
                      + theta_log_softmax[1].val_.val()
                      + theta_log_softmax[0].val_.val(),
                  stan::math::categorical_logit_log(ms, theta).val_.val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val() + theta_log_softmax[1].d_.val()
                      + theta_log_softmax[0].d_.val(),
                  stan::math::categorical_logit_log(ms, theta).d_.val());

  std::vector<vector_fv> theta_arr(3);
  theta_arr[0] = (vector_fv(3) << 2.6, 1.8, 6.5).finished();
  theta_arr[1] = (vector_fv(3) << 5.2, -1.0, 2.1).finished();
  theta_arr[2] = (vector_fv(3) << -2.3, 4.4, -3.8).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_ = i;
    theta_arr[1][i].d_ = i + i;
    theta_arr[2][i].d_ = i * 2;
  }

  std::vector<vector_fv> theta_arr_log_softmax = log_softmax(theta_arr);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].val_.val()
                      + theta_arr_log_softmax[1][1].val_.val()
                      + theta_arr_log_softmax[2][0].val_.val(),
                  stan::math::categorical_logit_log(ms, theta_arr).val_.val());
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].d_.val()
                      + theta_arr_log_softmax[1][1].d_.val()
                      + theta_arr_log_softmax[2][0].d_.val(),
                  stan::math::categorical_logit_log(ms, theta_arr).d_.val());

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].val_.val()
                      + theta_arr_log_softmax[1][1].val_.val()
                      + theta_arr_log_softmax[2][1].val_.val(),
                  stan::math::categorical_logit_log(2, theta_arr).val_.val());
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].d_.val()
                      + theta_arr_log_softmax[1][1].d_.val()
                      + theta_arr_log_softmax[2][1].d_.val(),
                  stan::math::categorical_logit_log(2, theta_arr).d_.val());
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_var) {
  using stan::math::vector_ffv;
  using stan::math::log_softmax;

  vector_ffv theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val() = i;
  vector_ffv theta_log_softmax = log_softmax(theta);
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val().val(),
                  stan::math::categorical_logit_log(1, theta).val_.val().val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].val_.val().val(),
                  stan::math::categorical_logit_log(2, theta).val_.val().val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].val_.val().val(),
                  stan::math::categorical_logit_log(3, theta).val_.val().val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val().val(),
                  stan::math::categorical_logit_log(1, theta).d_.val().val());
  EXPECT_FLOAT_EQ(theta_log_softmax[1].d_.val().val(),
                  stan::math::categorical_logit_log(2, theta).d_.val().val());
  EXPECT_FLOAT_EQ(theta_log_softmax[2].d_.val().val(),
                  stan::math::categorical_logit_log(3, theta).d_.val().val());
}

TEST(ProbDistributionsCategoricalLogit, fvar_fvar_var_vectorized) {
  using stan::math::vector_ffv;
  using stan::math::log_softmax;

  vector_ffv theta(3);
  theta << -1, 2, -10;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val() = i;

  std::vector<int> ns(0);
  EXPECT_FLOAT_EQ(0.0,
                  stan::math::categorical_logit_log(ns, theta).val_.val()
                                                                   .val());

  vector_ffv theta_log_softmax = log_softmax(theta);

  std::vector<int> ms(3);
  ms[0] = 1;
  ms[1] = 2;
  ms[2] = 1;
  EXPECT_FLOAT_EQ(theta_log_softmax[0].val_.val().val()
                      + theta_log_softmax[1].val_.val().val()
                      + theta_log_softmax[0].val_.val().val(),
                  stan::math::categorical_logit_log(ms, theta).val_.val()
                                                                   .val());
  EXPECT_FLOAT_EQ(theta_log_softmax[0].d_.val().val()
                      + theta_log_softmax[1].d_.val().val()
                      + theta_log_softmax[0].d_.val().val(),
                  stan::math::categorical_logit_log(ms, theta).d_.val().val());

  std::vector<vector_ffv> theta_arr(3);
  theta_arr[0] = (vector_ffv(3) << 1, 4, 3).finished();
  theta_arr[1] = (vector_ffv(3) << -1, -3, 2.1).finished();
  theta_arr[2] = (vector_ffv(3) << 0.3, 0.4, -2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_.val() = i;
    theta_arr[1][i].d_.val() = i * 3;
    theta_arr[2][i].d_.val() = i / 2;
  }

  std::vector<vector_ffv> theta_arr_log_softmax = log_softmax(theta_arr);

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].val_.val().val()
                      + theta_arr_log_softmax[1][1].val_.val().val()
                      + theta_arr_log_softmax[2][0].val_.val().val(),
                  stan::math::categorical_logit_log(ms, theta_arr).val_.val()
                                                                       .val());
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][0].d_.val().val()
                      + theta_arr_log_softmax[1][1].d_.val().val()
                      + theta_arr_log_softmax[2][0].d_.val().val(),
                  stan::math::categorical_logit_log(ms, theta_arr).d_.val()
                                                                     .val());

  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].val_.val().val()
                      + theta_arr_log_softmax[1][1].val_.val().val()
                      + theta_arr_log_softmax[2][1].val_.val().val(),
                  stan::math::categorical_logit_log(2, theta_arr).val_.val()
                                                                      .val());
  EXPECT_FLOAT_EQ(theta_arr_log_softmax[0][1].d_.val().val()
                      + theta_arr_log_softmax[1][1].d_.val().val()
                      + theta_arr_log_softmax[2][1].d_.val().val(),
                  stan::math::categorical_logit_log(2, theta_arr).d_.val()
                                                                    .val());
}
