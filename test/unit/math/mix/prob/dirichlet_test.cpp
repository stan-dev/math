#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

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
  EXPECT_FLOAT_EQ(0.99344212, stan::math::dirichlet_log(theta, alpha).d_.val());

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

TEST(ProbDistributions, fvar_varVectorised) {
  using stan::math::dirichlet_log;
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<var>, Dynamic, 1> theta1(3, 1), theta2(3, 1), theta3(3, 1);
  theta1 << 0.2, 0.3, 0.5;
  theta2 << 0.1, 0.8, 0.1;
  theta3 << 0.6, 0.1, 0.3;

  Matrix<fvar<var>, Dynamic, 1> alpha1(3, 1), alpha2(3, 1), alpha3(3, 1);
  alpha1 << 1.0, 1.0, 1.0;
  alpha2 << 6.2, 3.5, 9.1;
  alpha3 << 2.5, 7.4, 6.1;

  std::vector<Matrix<fvar<var>, Dynamic, 1>> theta_vec(3);
  theta_vec[0] = theta1;
  theta_vec[1] = theta2;
  theta_vec[2] = theta3;

  std::vector<Matrix<fvar<var>, Dynamic, 1>> alpha_vec(3);
  alpha_vec[0] = alpha1;
  alpha_vec[1] = alpha2;
  alpha_vec[2] = alpha3;

  Matrix<fvar<var>, Dynamic, 1> result(3);
  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta2, alpha2);
  result[2] = dirichlet_log(theta3, alpha3);

  fvar<var> out = dirichlet_log(theta_vec, alpha_vec);

  EXPECT_FLOAT_EQ(result.val().val().sum(), out.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().sum(), out.d_.val());

  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta2, alpha1);
  result[2] = dirichlet_log(theta3, alpha1);

  out = dirichlet_log(theta_vec, alpha1);

  EXPECT_FLOAT_EQ(result.val().val().sum(), out.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().sum(), out.d_.val());

  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta1, alpha2);
  result[2] = dirichlet_log(theta1, alpha3);

  out = dirichlet_log(theta1, alpha_vec);

  EXPECT_FLOAT_EQ(result.val().val().sum(), out.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().sum(), out.d_.val());
}

TEST(ProbDistributions, fvar_fvar_var) {
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<fvar<var>>, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  Matrix<fvar<fvar<var>>, Dynamic, 1> alpha(3, 1);
  alpha << 1.0, 1.0, 1.0;
  for (int i = 0; i < 3; i++) {
    theta(i).d_ = 1.0;
    alpha(i).d_ = 1.0;
  }

  EXPECT_FLOAT_EQ(0.6931472,
                  stan::math::dirichlet_log(theta, alpha).val_.val_.val());
  EXPECT_FLOAT_EQ(0.99344212,
                  stan::math::dirichlet_log(theta, alpha).d_.val_.val());

  Matrix<fvar<fvar<var>>, Dynamic, 1> theta2(4, 1);
  theta2 << 0.01, 0.01, 0.8, 0.18;
  Matrix<fvar<fvar<var>>, Dynamic, 1> alpha2(4, 1);
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

TEST(ProbDistributions, fvar_fvar_varVectorised) {
  using stan::math::dirichlet_log;
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<fvar<var>>, Dynamic, 1> theta1(3, 1), theta2(3, 1), theta3(3, 1);
  theta1 << 0.2, 0.3, 0.5;
  theta2 << 0.1, 0.8, 0.1;
  theta3 << 0.6, 0.1, 0.3;

  Matrix<fvar<fvar<var>>, Dynamic, 1> alpha1(3, 1), alpha2(3, 1), alpha3(3, 1);
  alpha1 << 1.0, 1.0, 1.0;
  alpha2 << 6.2, 3.5, 9.1;
  alpha3 << 2.5, 7.4, 6.1;

  std::vector<Matrix<fvar<fvar<var>>, Dynamic, 1>> theta_vec(3);
  theta_vec[0] = theta1;
  theta_vec[1] = theta2;
  theta_vec[2] = theta3;

  std::vector<Matrix<fvar<fvar<var>>, Dynamic, 1>> alpha_vec(3);
  alpha_vec[0] = alpha1;
  alpha_vec[1] = alpha2;
  alpha_vec[2] = alpha3;

  Matrix<fvar<fvar<var>>, Dynamic, 1> result(3);
  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta2, alpha2);
  result[2] = dirichlet_log(theta3, alpha3);

  fvar<fvar<var>> out = dirichlet_log(theta_vec, alpha_vec);

  EXPECT_FLOAT_EQ(result.val().val().val().sum(), out.val_.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().val().sum(), out.d_.val_.val());

  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta2, alpha1);
  result[2] = dirichlet_log(theta3, alpha1);

  out = dirichlet_log(theta_vec, alpha1);

  EXPECT_FLOAT_EQ(result.val().val().val().sum(), out.val_.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().val().sum(), out.d_.val_.val());

  result[0] = dirichlet_log(theta1, alpha1);
  result[1] = dirichlet_log(theta1, alpha2);
  result[2] = dirichlet_log(theta1, alpha3);

  out = dirichlet_log(theta1, alpha_vec);

  EXPECT_FLOAT_EQ(result.val().val().val().sum(), out.val_.val_.val());
  EXPECT_FLOAT_EQ(result.d().val().val().sum(), out.d_.val_.val());
}
