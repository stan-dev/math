#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

TEST(ProbDistributionsCategorical, fvar_var) {
  using stan::math::vector_fv;

  vector_fv theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  EXPECT_FLOAT_EQ(std::log(0.3),
                  stan::math::categorical_log(1, theta).val_.val());
  EXPECT_FLOAT_EQ(std::log(0.5),
                  stan::math::categorical_log(2, theta).val_.val());
  EXPECT_FLOAT_EQ(std::log(0.2),
                  stan::math::categorical_log(3, theta).val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.3, stan::math::categorical_log(1, theta).d_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.5, stan::math::categorical_log(2, theta).d_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.2, stan::math::categorical_log(3, theta).d_.val());
}
TEST(ProbDistributionsCategorical, fvar_var_vector) {
  using stan::math::vector_fv;

  vector_fv theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_log(xs, theta).val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.3 + 1.0 / 0.2 + 1.0 / 0.3,
                  stan::math::categorical_log(xs, theta).d_.val());

  std::vector<vector_fv> theta_arr(3);
  theta_arr[0] = (vector_fv(3) << 0.3, 0.4, 0.3).finished();
  theta_arr[1] = (vector_fv(3) << 0.1, 0.3, 0.6).finished();
  theta_arr[2] = (vector_fv(3) << 0.6, 0.2, 0.2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_ = i;
    theta_arr[1][i].d_ = i * 3;
    theta_arr[2][i].d_ = i / 2;
  }

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val()
                      + log(theta_arr[1][2]).val_.val()
                      + log(theta_arr[2][0]).val_.val(),
                  stan::math::categorical_log(xs, theta_arr).val_.val());
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val() + log(theta_arr[1][2]).d_.val()
                      + log(theta_arr[2][0]).d_.val(),
                  stan::math::categorical_log(xs, theta_arr).d_.val());

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val()
                      + log(theta_arr[1][0]).val_.val()
                      + log(theta_arr[2][0]).val_.val(),
                  stan::math::categorical_log(1, theta_arr).val_.val());
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val() + log(theta_arr[1][0]).d_.val()
                      + log(theta_arr[2][0]).d_.val(),
                  stan::math::categorical_log(1, theta_arr).d_.val());
}

TEST(ProbDistributionsCategorical, fvar_fvar_var) {
  using stan::math::vector_ffv;

  vector_ffv theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  EXPECT_FLOAT_EQ(std::log(0.3),
                  stan::math::categorical_log(1, theta).val_.val_.val());
  EXPECT_FLOAT_EQ(std::log(0.5),
                  stan::math::categorical_log(2, theta).val_.val_.val());
  EXPECT_FLOAT_EQ(std::log(0.2),
                  stan::math::categorical_log(3, theta).val_.val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.3,
                  stan::math::categorical_log(1, theta).d_.val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.5,
                  stan::math::categorical_log(2, theta).d_.val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.2,
                  stan::math::categorical_log(3, theta).d_.val_.val());
}
TEST(ProbDistributionsCategorical, fvar_fvar_var_vector) {
  using stan::math::vector_ffv;

  vector_ffv theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_log(xs, theta).val_.val_.val());
  EXPECT_FLOAT_EQ(1.0 / 0.3 + 1.0 / 0.2 + 1.0 / 0.3,
                  stan::math::categorical_log(xs, theta).d_.val_.val());

  std::vector<vector_ffv> theta_arr(3);
  theta_arr[0] = (vector_ffv(3) << 0.8, 0.1, 0.1).finished();
  theta_arr[1] = (vector_ffv(3) << 0.2, 0.6, 0.2).finished();
  theta_arr[2] = (vector_ffv(3) << 0.2, 0.1, 0.7).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_.val_ = i;
    theta_arr[1][i].d_.val_ = i * 3;
    theta_arr[2][i].d_.val_ = i / 2;
  }

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val_.val()
                      + log(theta_arr[1][2]).val_.val_.val()
                      + log(theta_arr[2][0]).val_.val_.val(),
                  stan::math::categorical_log(xs, theta_arr).val_.val_.val());
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val_.val()
                      + log(theta_arr[1][2]).d_.val_.val()
                      + log(theta_arr[2][0]).d_.val_.val(),
                  stan::math::categorical_log(xs, theta_arr).d_.val_.val());

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val_.val()
                      + log(theta_arr[1][0]).val_.val_.val()
                      + log(theta_arr[2][0]).val_.val_.val(),
                  stan::math::categorical_log(1, theta_arr).val_.val_.val());
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val_.val()
                      + log(theta_arr[1][0]).d_.val_.val()
                      + log(theta_arr[2][0]).d_.val_.val(),
                  stan::math::categorical_log(1, theta_arr).d_.val_.val());
}
