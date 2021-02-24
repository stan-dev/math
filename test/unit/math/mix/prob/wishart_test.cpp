#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsWishart, matvar) {
  auto f = [](const auto& y, const auto& dof, const auto& sigma) {
    auto y_sym = stan::math::multiply(0.5, y + y.transpose());
    auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
    return stan::math::wishart_lpdf(y_sym, dof, sigma_sym);
  };

  double dof = 3.4;
  Eigen::MatrixXd y11(1, 1);
  y11 << 1;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  stan::test::expect_ad(f, y11, dof, Sigma11);
  stan::test::expect_ad_matvar(f, y11, dof, Sigma11);

  Eigen::MatrixXd y00(0, 0);
  Eigen::MatrixXd Sigma00(0, 0);
  stan::test::expect_ad(f, y00, dof, Sigma00);
  stan::test::expect_ad_matvar(f, y00, dof, Sigma00);

  Eigen::MatrixXd y22(2, 2);
  y22 << 1.0, 0.1, 0.1, 2.0;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 2.0, 0.5, 0.5, 1.1;
  stan::test::expect_ad(f, y22, dof, Sigma22);
  stan::test::expect_ad_matvar(f, y22, dof, Sigma22);

  // Error sizes
  stan::test::expect_ad(f, y00, dof, Sigma11);
  stan::test::expect_ad(f, y11, dof, Sigma00);
  stan::test::expect_ad_matvar(f, y00, dof, Sigma11);
  stan::test::expect_ad_matvar(f, y11, dof, Sigma00);
}

TEST(ProbDistributionsWishart, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<var>, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.20661, -11.20661, 112.94139;

  for (int i = 0; i < 4; i++) {
    Sigma(i).d_ = 1.0;
    Y(i).d_ = 1.0;
  }

  unsigned int dof = 3;
  // computed with MCMCpack in R
  double lp = log(8.658e-07);

  EXPECT_NEAR(lp, stan::math::wishart_log(Y, dof, Sigma).val_.val(), 0.01);
  EXPECT_NEAR(-0.76893887, stan::math::wishart_log(Y, dof, Sigma).d_.val(),
              0.01);

  stan::math::recover_memory();
}

TEST(ProbDistributionsWishart, fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(2, 2);
  Sigma << 1.848220, 1.899623, 1.899623, 12.751941;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Y(2, 2);
  Y << 2.011108, -11.206611, -11.206611, 112.94139;

  for (int i = 0; i < 4; i++) {
    Sigma(i).d_.val_ = 1.0;
    Y(i).d_.val_ = 1.0;
  }

  unsigned int dof = 3;
  // computed with MCMCpack in R
  double lp = log(8.658e-07);

  EXPECT_NEAR(lp, stan::math::wishart_log(Y, dof, Sigma).val_.val_.val(), 0.01);
  EXPECT_NEAR(-0.76893887, stan::math::wishart_log(Y, dof, Sigma).d_.val_.val(),
              0.01);

  stan::math::recover_memory();
}
