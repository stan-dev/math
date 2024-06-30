#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, ProbDistributionsMatrixNormal_matvar) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma, const auto& D) {
    auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
    auto D_sym = stan::math::multiply(0.5, D + D.transpose());
    return stan::math::matrix_normal_prec_lpdf(y, mu, sigma_sym, D_sym);
  };

  auto f_const_y = [](const auto& y) {
    return [&y](const auto& mu, const auto& sigma, const auto& D) {
      auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
      auto D_sym = stan::math::multiply(0.5, D + D.transpose());
      return stan::math::matrix_normal_prec_lpdf(y, mu, sigma_sym, D_sym);
    };
  };

  auto f_const_D = [](const auto& D) {
    return [&D](const auto& y, const auto& mu, const auto& sigma) {
      auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
      auto D_sym = stan::math::multiply(0.5, D + D.transpose());
      return stan::math::matrix_normal_prec_lpdf(y, mu, sigma_sym, D_sym);
    };
  };

  Eigen::VectorXd y1(1);
  y1 << 1;
  Eigen::VectorXd mu1(1);
  mu1 << 3.4;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  Eigen::MatrixXd D11(1, 1);
  D11 << 1.1;
  stan::test::expect_ad(f_const_y(y1), mu1, Sigma11, D11);
  stan::test::expect_ad(f_const_D(D11), y1, mu1, Sigma11);
  stan::test::expect_ad_matvar(f, y1, mu1, Sigma11, D11);

  Eigen::VectorXd y0(0);
  Eigen::VectorXd mu0(0);
  Eigen::MatrixXd Sigma00(0, 0);
  Eigen::MatrixXd D00(0, 0);
  stan::test::expect_ad(f_const_y(y0), mu0, Sigma00, D00);
  stan::test::expect_ad(f_const_D(D00), y0, mu0, Sigma00);
  stan::test::expect_ad_matvar(f, y0, mu0, Sigma00, D00);

  Eigen::VectorXd y2(2);
  y2 << 1.0, 0.1;
  Eigen::VectorXd mu2(2);
  mu2 << 0.1, 2.0;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 2.0, 0.5, 0.5, 1.1;
  Eigen::MatrixXd D22(2, 2);
  Sigma22 << 1.1, 0.5, 0.5, 1.2;
  stan::test::expect_ad(f_const_y(y2), mu2, Sigma22, D22);
  stan::test::expect_ad(f_const_D(D22), y2, mu2, Sigma22);
  stan::test::expect_ad_matvar(f, y2, mu2, Sigma22, D22);

  // Error sizes
  stan::test::expect_ad(f_const_y(y0), mu0, Sigma11, D11);
  stan::test::expect_ad(f_const_D(D11), y0, mu0, Sigma11);
  stan::test::expect_ad(f_const_y(y1), mu0, Sigma00, D11);
  stan::test::expect_ad(f_const_D(D11), y1, mu0, Sigma00);
  stan::test::expect_ad_matvar(f, y0, mu0, Sigma11, D11);
  stan::test::expect_ad_matvar(f, y1, mu1, Sigma00, D11);
}

TEST_F(AgradRev, ProbDistributionsMatrixNormal_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<var>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<var>, Dynamic, Dynamic> mu(3, 5);
  mu.setZero();

  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  Matrix<fvar<var>, Dynamic, Dynamic> D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++) {
      D(i, j).d_ = 1.0;
      if (i < 3) {
        mu(i, j).d_ = 1.0;
        y(i, j).d_ = 1.0;
        if (j < 3)
          Sigma(i, j).d_ = 1.0;
      }
    }

  fvar<var> lp_ref = stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D);
  EXPECT_FLOAT_EQ(-2132.07482, lp_ref.val_.val());
  EXPECT_FLOAT_EQ(-2075.1274, lp_ref.d_.val());

  stan::math::recover_memory();
}

TEST_F(AgradRev, ProbDistributionsMatrixNormal_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;

  Matrix<fvar<fvar<var>>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<fvar<var>>, Dynamic, Dynamic> mu(3, 5);
  mu.setZero();

  Matrix<fvar<fvar<var>>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 1.0, 0.5, 0.1, 0.5, 1.0, 0.2, 0.1, 0.2, 1.0;

  Matrix<fvar<fvar<var>>, Dynamic, Dynamic> D(5, 5);
  D << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 1.0,
      0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++) {
      D(i, j).d_.val_ = 1.0;
      if (i < 3) {
        mu(i, j).d_.val_ = 1.0;
        y(i, j).d_.val_ = 1.0;
        if (j < 3)
          Sigma(i, j).d_.val_ = 1.0;
      }
    }

  fvar<fvar<var>> lp_ref = stan::math::matrix_normal_prec_lpdf(y, mu, Sigma, D);
  EXPECT_FLOAT_EQ(-2132.07482, lp_ref.val_.val_.val());
  EXPECT_FLOAT_EQ(-2075.1274, lp_ref.d_.val_.val());

  stan::math::recover_memory();
}
