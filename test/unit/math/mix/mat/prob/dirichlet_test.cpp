#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;


TEST(ProbDistributions, fvar_var) {
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<var>, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  Matrix<fvar<var>, Dynamic, 1> alpha(3, 1);
  alpha << 1.0, 1.0, 1.0;
  for (int i = 0; i < 3; i++) {
    theta(i).d_ = 1.0;
    alpha(i).d_ = 1.0;
  }

  EXPECT_FLOAT_EQ(0.6931472,
                  stan::math::dirichlet_log(theta, alpha).val_.val());
  EXPECT_FLOAT_EQ(0.99344212,
                  stan::math::dirichlet_log(theta, alpha).d_.val());

  Matrix<fvar<var>, Dynamic, 1> theta2(4, 1);
  theta2 << 0.01, 0.01, 0.8, 0.18;
  Matrix<fvar<var>, Dynamic, 1> alpha2(4, 1);
  alpha2 << 10.5, 11.5, 19.3, 5.1;
  for (int i = 0; i < 3; i++) {
    theta2(i).d_ = 1.0;
    alpha2(i).d_ = 1.0;
  }

  EXPECT_FLOAT_EQ(-43.40045,
                  stan::math::dirichlet_log(theta2, alpha2).val_.val());
  EXPECT_FLOAT_EQ(2017.2858,
                  stan::math::dirichlet_log(theta2, alpha2).d_.val());
}

TEST(ProbDistributions, fvar_fvar_var) {
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<fvar<var> >, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  Matrix<fvar<fvar<var> >, Dynamic, 1> alpha(3, 1);
  alpha << 1.0, 1.0, 1.0;
  for (int i = 0; i < 3; i++) {
    theta(i).d_ = 1.0;
    alpha(i).d_ = 1.0;
  }

  EXPECT_FLOAT_EQ(0.6931472,
                  stan::math::dirichlet_log(theta, alpha).val_.val_.val());
  EXPECT_FLOAT_EQ(0.99344212,
                  stan::math::dirichlet_log(theta, alpha).d_.val_.val());

  Matrix<fvar<fvar<var> >, Dynamic, 1> theta2(4, 1);
  theta2 << 0.01, 0.01, 0.8, 0.18;
  Matrix<fvar<fvar<var> >, Dynamic, 1> alpha2(4, 1);
  alpha2 << 10.5, 11.5, 19.3, 5.1;
  for (int i = 0; i < 3; i++) {
    theta2(i).d_ = 1.0;
    alpha2(i).d_ = 1.0;
  }

  EXPECT_FLOAT_EQ(-43.40045,
                  stan::math::dirichlet_log(theta2, alpha2).val_.val_.val());
  EXPECT_FLOAT_EQ(2017.2858,
                  stan::math::dirichlet_log(theta2, alpha2).d_.val_.val());
}
