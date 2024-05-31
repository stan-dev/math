#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsWishartCholesky, matvar) {
  auto f = [](const auto& Y, const auto& dof, const auto& Sigma) {
    auto symmetric_Y = ((Y + Y.transpose()) * 0.5).eval();
    auto symmetric_Sigma = ((Sigma + Sigma.transpose()) * 0.5).eval();

    auto L_Y = stan::math::cholesky_decompose(symmetric_Y);
    auto L_S = stan::math::cholesky_decompose(symmetric_Sigma);

    return stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S);
  };

  double dof = 3.4;
  Eigen::MatrixXd Y11(1, 1);
  Y11 << 1;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  stan::test::expect_ad(f, Y11, dof, Sigma11);
  stan::test::expect_ad_matvar(f, Y11, dof, Sigma11);

  Eigen::MatrixXd Y00(0, 0);
  Eigen::MatrixXd Sigma00(0, 0);
  stan::test::expect_ad(f, Y00, dof, Sigma00);
  stan::test::expect_ad_matvar(f, Y00, dof, Sigma00);

  Eigen::MatrixXd Y22(2, 2);
  Y22 << 1.0, 0.1, 0.1, 1.5;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 1.5, 0.5, 0.5, 1.1;
  stan::test::expect_ad(f, Y22, dof, Sigma22);
  stan::test::expect_ad_matvar(f, Y22, dof, Sigma22);

  // Error sizes
  stan::test::expect_ad(f, Y00, dof, Sigma11);
  stan::test::expect_ad(f, Y11, dof, Sigma00);
  stan::test::expect_ad_matvar(f, Y00, dof, Sigma11);
  stan::test::expect_ad_matvar(f, Y11, dof, Sigma00);
}

TEST(ProbDistributionsWishartCholesky, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<var>, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.20661, -11.20661, 112.94139;

  auto L_Y = stan::math::cholesky_decompose(Y);
  auto L_S = stan::math::cholesky_decompose(Sigma);

  for (int i = 0; i < 4; i++) {
    L_Y(i).d_ = 1.0;
    L_S(i).d_ = 1.0;
  }

  unsigned int dof = 3;
  // log absolute determinant of the change of variables from Y -> LL'
  // see Theorem 2.1.9 in Muirhead, Aspects of Multivariate Statistical Theory
  // computed with MCMCpack in R
  double lp = -13.95961162478571360 + 4.04590909757044237;

  EXPECT_NEAR(lp, stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S).val_.val(),
              1e-9);
  EXPECT_NEAR(2.3751897169452936,
              stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S).d_.val(), 1e-9);

  stan::math::recover_memory();
}

TEST(ProbDistributionsWishartCholesky, fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  auto L_Y = stan::math::cholesky_decompose(Y);
  auto L_S = stan::math::cholesky_decompose(Sigma);

  for (int i = 0; i < 4; i++) {
    L_Y(i).d_ = 1.0;
    L_S(i).d_ = 1.0;
  }
  unsigned int dof = 3;
  // log absolute determinant of the change of variables from Y -> LL'
  // see Theorem 2.1.9 in Muirhead, Aspects of Multivariate Statistical Theory
  // computed with MCMCpack in R
  double lp = -13.95961162478571360 + 4.04590909757044237;

  EXPECT_NEAR(lp,
              stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S).val_.val_.val(),
              1e-6);
  EXPECT_NEAR(2.3751897169452936,
              stan::math::wishart_cholesky_lpdf(L_Y, dof, L_S).d_.val_.val(),
              1e-6);

  stan::math::recover_memory();
}
