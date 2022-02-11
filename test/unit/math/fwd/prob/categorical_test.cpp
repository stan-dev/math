#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

TEST(ProbDistributionsCategorical, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<double>, Dynamic, 1> theta(3, 1);
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
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<double>, Dynamic, 1> theta(3, 1);
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
}

TEST(ProbDistributionsCategorical, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta(3, 1);
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
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta(3, 1);
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
}
