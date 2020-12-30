#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

TEST(ProbDistributionsCategorical, fvar_double) {
  using stan::math::vector_fd;
  vector_fd theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  EXPECT_FLOAT_EQ(std::log(0.3), stan::math::categorical_log(1, theta).val_);
  EXPECT_FLOAT_EQ(std::log(0.5), stan::math::categorical_log(2, theta).val_);
  EXPECT_FLOAT_EQ(std::log(0.2), stan::math::categorical_log(3, theta).val_);
  EXPECT_FLOAT_EQ(1.0 / 0.3, stan::math::categorical_log(1, theta).d_);
  EXPECT_FLOAT_EQ(1.0 / 0.5, stan::math::categorical_log(2, theta).d_);
  EXPECT_FLOAT_EQ(1.0 / 0.2, stan::math::categorical_log(3, theta).d_);
}
TEST(ProbDistributionsCategorical, fvar_double_vector) {
  using stan::math::vector_fd;
  vector_fd theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_log(xs, theta).val_);
  EXPECT_FLOAT_EQ(1.0 / 0.3 + 1.0 / 0.2 + 1.0 / 0.3,
                  stan::math::categorical_log(xs, theta).d_);

  std::vector<vector_fd> theta_arr(3);
  theta_arr[0] = (vector_fd(3) << 0.3, 0.4, 0.3).finished();
  theta_arr[1] = (vector_fd(3) << 0.1, 0.3, 0.6).finished();
  theta_arr[2] = (vector_fd(3) << 0.6, 0.2, 0.2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_ = i;
    theta_arr[1][i].d_ = i * 3;
    theta_arr[2][i].d_ = i / 2;
  }

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_ + log(theta_arr[1][2]).val_
                  + log(theta_arr[2][0]).val_,
                  stan::math::categorical_log(xs, theta_arr).val_);
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_ + log(theta_arr[1][2]).d_
                  + log(theta_arr[2][0]).d_,
                  stan::math::categorical_log(xs, theta_arr).d_);

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_ + log(theta_arr[1][0]).val_
                  + log(theta_arr[2][0]).val_,
                  stan::math::categorical_log(1, theta_arr).val_);
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_ + log(theta_arr[1][0]).d_
                  + log(theta_arr[2][0]).d_,
                  stan::math::categorical_log(1, theta_arr).d_);
}

TEST(ProbDistributionsCategorical, fvar_fvar_double) {
  using stan::math::vector_ffd;
  vector_ffd theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  EXPECT_FLOAT_EQ(std::log(0.3),
                  stan::math::categorical_log(1, theta).val_.val_);
  EXPECT_FLOAT_EQ(std::log(0.5),
                  stan::math::categorical_log(2, theta).val_.val_);
  EXPECT_FLOAT_EQ(std::log(0.2),
                  stan::math::categorical_log(3, theta).val_.val_);
  EXPECT_FLOAT_EQ(1.0 / 0.3, stan::math::categorical_log(1, theta).d_.val_);
  EXPECT_FLOAT_EQ(1.0 / 0.5, stan::math::categorical_log(2, theta).d_.val_);
  EXPECT_FLOAT_EQ(1.0 / 0.2, stan::math::categorical_log(3, theta).d_.val_);
}
TEST(ProbDistributionsCategorical, fvar_fvar_double_vector) {
  using stan::math::vector_ffd;
  vector_ffd theta(3);
  theta << 0.3, 0.5, 0.2;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  std::vector<int> xs(3);
  xs[0] = 1;
  xs[1] = 3;
  xs[2] = 1;

  EXPECT_FLOAT_EQ(log(0.3) + log(0.2) + log(0.3),
                  stan::math::categorical_log(xs, theta).val_.val_);
  EXPECT_FLOAT_EQ(1.0 / 0.3 + 1.0 / 0.2 + 1.0 / 0.3,
                  stan::math::categorical_log(xs, theta).d_.val_);

  std::vector<vector_ffd> theta_arr(3);
  theta_arr[0] = (vector_ffd(3) << 0.3, 0.4, 0.3).finished();
  theta_arr[1] = (vector_ffd(3) << 0.1, 0.3, 0.6).finished();
  theta_arr[2] = (vector_ffd(3) << 0.6, 0.2, 0.2).finished();

  for (int i = 0; i < 3; i++) {
    theta_arr[0][i].d_.val_ = i;
    theta_arr[1][i].d_.val_ = i * 3;
    theta_arr[2][i].d_.val_ = i / 2;
  }

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val_
                    + log(theta_arr[1][2]).val_.val_
                    + log(theta_arr[2][0]).val_.val_,
                  stan::math::categorical_log(xs, theta_arr).val_.val_);
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val_
                    + log(theta_arr[1][2]).d_.val_
                    + log(theta_arr[2][0]).d_.val_,
                  stan::math::categorical_log(xs, theta_arr).d_.val_);

  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).val_.val_
                    + log(theta_arr[1][0]).val_.val_
                    + log(theta_arr[2][0]).val_.val_,
                  stan::math::categorical_log(1, theta_arr).val_.val_);
  EXPECT_FLOAT_EQ(log(theta_arr[0][0]).d_.val_
                    + log(theta_arr[1][0]).d_.val_
                    + log(theta_arr[2][0]).d_.val_,
                  stan::math::categorical_log(1, theta_arr).d_.val_);
}
