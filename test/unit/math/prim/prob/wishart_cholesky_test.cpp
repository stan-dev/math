#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions.hpp>

TEST(ProbDistributionsWishartCholesky, wishart_symmetry) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::wishart_cholesky_lpdf;

  MatrixXd Sigma(4, 4);
  MatrixXd Y(4, 4);

  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570,
      49.215769, -15.454186, -14.474830, 49.215769, 60.08987, -20.48108,
      4.395895, -15.454186, -20.48108, 7.885833;

  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.586685,
      0.1055911, -3.1129955, -3.586685, 1.4482736;

  unsigned int dof = 5;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  EXPECT_NO_THROW(wishart_cholesky_lpdf(L_Y, dof, L_S));
}

TEST(ProbDistributionsWishartCholesky, wishart_pos_def) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::wishart_cholesky_lpdf;

  MatrixXd Sigma(2, 2);
  MatrixXd Sigma_non_pos_def(2, 2);
  MatrixXd Y(2, 2);
  MatrixXd Y_non_pos_def(2, 2);

  Sigma << 1, 0, 0, 1;
  Sigma_non_pos_def << -1, 0, 0, 1;
  Y << 1, 0, 0, 1;
  Y_non_pos_def << -1, 0, 0, 1;

  unsigned int dof = 5;

  EXPECT_NO_THROW(wishart_cholesky_lpdf(Y, dof, Sigma));
  EXPECT_THROW(wishart_cholesky_lpdf(Y_non_pos_def, dof, Sigma),
               std::domain_error);
  EXPECT_THROW(wishart_cholesky_lpdf(Y, dof, Sigma_non_pos_def),
               std::domain_error);
}

TEST(ProbDistributionsWishartCholesky, 0x0) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Sigma(0, 0);
  MatrixXd Y(0, 0);

  unsigned int dof = 3;
  EXPECT_THROW(stan::math::wishart_cholesky_lpdf(Y, dof, Sigma),
               std::domain_error);

  unsigned int dof0 = 0;
  EXPECT_THROW(stan::math::wishart_cholesky_lpdf(Y, dof0, Sigma),
               std::domain_error);
}

TEST(ProbDistributionsWishartCholesky, dof_0) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  MatrixXd Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  unsigned int dof = 0;
  EXPECT_THROW(stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S),
               std::domain_error);
}

TEST(ProbDistributionsWishartCholesky, 1x1) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Sigma(1, 1);
  Sigma << 1;

  MatrixXd Y(1, 1);
  Y << 2.011108;

  double dof = 0.1;

  EXPECT_NO_THROW(stan::math::wishart_cholesky_lpdf(Y, dof, Sigma));
}

TEST(ProbDistributionsWishartCholesky, 2x2) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  MatrixXd Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  unsigned int dof = 3;
  // log absolute determinant of the change of variables from Y -> LL'
  // see Theorem 2.1.9 in Muirhead, Aspects of Multivariate Statistical Theory
  double log_jac = 2 * stan::math::LOG_TWO;

  for (int i = 0; i < 2; i++) {
    log_jac += (2 - i) * log(L_Y(i, i));
  }
  // computed with MCMCpack in R
  double lp = -13.9596117200 + log_jac;

  EXPECT_NEAR(lp, stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S), 1e-9);
}

TEST(ProbDistributionsWishartCholesky, 4x4) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Y(4, 4);
  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570, 49.21577,
      -18.454186, -14.474830, 49.21577, 60.08987, -21.48108, 4.395895,
      -18.454186, -21.48108, 7.885833;

  MatrixXd Sigma(4, 4);
  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.5866852, 1.4482736;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  unsigned int dof = 4;
  double log_jac = 4 * stan::math::LOG_TWO;
  unsigned int n = 4;
  for (int i = 0; i < n; i++) {
    log_jac += (n - i) * log(L_Y(i, i));
  }

  double log_p = -20.9420668184 + log_jac;
  EXPECT_NEAR(log_p, stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S), 1e-9);

  dof = 5;

  log_p = -22.6084823429 + log_jac;
  EXPECT_NEAR(log_p, stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S), 1e-9);
}

TEST(ProbDistributionsWishartCholesky, 2x2Propto) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  MatrixXd Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  unsigned int dof = 3;

  EXPECT_FLOAT_EQ(0.0, stan::math::wishart_cholesky_lpdf<true>(L_Y, dof, L_S));
}

TEST(ProbDistributionsWishartCholesky, 4x4Propto) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  MatrixXd Y(4, 4);
  Y << 7.988168, -9.555605, -14.47483, 4.395895, -9.555605, 44.750570,
      49.215769, -18.454186, -14.474830, 49.215769, 60.08987, -21.48108,
      4.395895, -18.454186, -21.48108, 7.885833;

  MatrixXd Sigma(4, 4);
  Sigma << 2.9983662, 0.2898776, -2.650523, 0.1055911, 0.2898776, 11.4803610,
      7.1579931, -3.1129955, -2.650523, 7.1579931, 11.676181, -3.5866852,
      0.1055911, -3.1129955, -3.5866852, 1.4482736;

  MatrixXd L_Y = Y.llt().matrixL();
  MatrixXd L_S = Sigma.llt().matrixL();

  double dof = 4;
  EXPECT_FLOAT_EQ(0, stan::math::wishart_cholesky_lpdf<true>(L_Y, dof, L_S));

  dof = 5;
  EXPECT_FLOAT_EQ(0, stan::math::wishart_cholesky_lpdf<true>(L_Y, dof, L_S));
}

TEST(ProbDistributionsWishartCholesky, error) {
  using Eigen::Dynamic;
  using Eigen::MatrixXd;
  using stan::math::wishart_cholesky_lpdf;
  double nu;

  nu = 1;
  EXPECT_NO_THROW(wishart_cholesky_lpdf(MatrixXd::Identity(1, 1), nu,
                                        MatrixXd::Identity(1, 1)));

  nu = 5;
  MatrixXd Sigma(2, 1);
  EXPECT_THROW(wishart_cholesky_lpdf(MatrixXd::Identity(1, 1), nu, Sigma),
               std::invalid_argument);

  nu = 5;
  MatrixXd Y(2, 1);
  EXPECT_THROW(wishart_cholesky_lpdf(Y, nu, MatrixXd::Identity(2, 2)),
               std::invalid_argument);

  nu = 5;
  EXPECT_THROW(wishart_cholesky_lpdf(MatrixXd::Identity(3, 3), nu,
                                     MatrixXd::Identity(2, 2)),
               std::invalid_argument);

  nu = 3;
  EXPECT_NO_THROW(wishart_cholesky_lpdf(MatrixXd::Identity(3, 3), nu,
                                        MatrixXd::Identity(3, 3)));
  nu = 2;
  EXPECT_THROW(wishart_cholesky_lpdf(MatrixXd::Identity(3, 3), nu,
                                     MatrixXd::Identity(3, 3)),
               std::domain_error);
}
