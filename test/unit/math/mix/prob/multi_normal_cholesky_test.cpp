#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, ProbDistributionsMultiNormalCholesky_matvar) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
    auto L = stan::math::cholesky_decompose(sigma_sym);
    return stan::math::multi_normal_cholesky_lpdf(y, mu, L);
  };

  Eigen::VectorXd y1(1);
  y1 << 1;
  Eigen::VectorXd mu1(1);
  mu1 << 3.4;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  stan::test::expect_ad(f, y1, mu1, Sigma11);

  Eigen::VectorXd y0(0);
  Eigen::VectorXd mu0(0);
  Eigen::MatrixXd Sigma00(0, 0);
  stan::test::expect_ad(f, y0, mu0, Sigma00);

  Eigen::VectorXd y2(2);
  y2 << 1.0, 0.1;
  Eigen::VectorXd mu2(2);
  mu2 << 0.1, 2.0;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 2.0, 0.5, 0.5, 1.1;
  stan::test::expect_ad(f, y2, mu2, Sigma22);

  Eigen::VectorXd y22(2);
  y22 << 0.4, 0.3;
  Eigen::VectorXd mu22(2);
  mu22 << 2.1, 1.0;
  std::vector<Eigen::VectorXd> y2s = {y2, y22};
  std::vector<Eigen::VectorXd> mu2s = {mu2, mu22};
  stan::test::expect_ad(f, y2, mu2s, Sigma22);
  stan::test::expect_ad(f, y2s, mu2, Sigma22);
  stan::test::expect_ad(f, y2s, mu2s, Sigma22);

  // Error sizes
  stan::test::expect_ad(f, y0, mu0, Sigma11);
  stan::test::expect_ad(f, y1, mu1, Sigma00);
  stan::test::expect_ad_matvar(f, y0, mu0, Sigma11);
  stan::test::expect_ad_matvar(f, y1, mu1, Sigma00);
}

TEST_F(AgradRev, ProbDistributionsMultiNormalCholesky_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  using std::vector;
  Matrix<fvar<var>, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<var>, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;

  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_ = 1.0;
  }

  Matrix<fvar<var>, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  EXPECT_FLOAT_EQ(-11.73908,
                  stan::math::multi_normal_cholesky_lpdf(y, mu, L).val_.val());
  EXPECT_FLOAT_EQ(0.54899865,
                  stan::math::multi_normal_cholesky_lpdf(y, mu, L).d_.val());
}

TEST_F(AgradRev, ProbDistributionsMultiNormalCholesky_fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  using std::vector;
  Matrix<fvar<fvar<var> >, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<fvar<var> >, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;

  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_.val_ = 1.0;
    mu(i).d_.val_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_.val_ = 1.0;
  }

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  EXPECT_FLOAT_EQ(
      -11.73908,
      stan::math::multi_normal_cholesky_lpdf(y, mu, L).val_.val_.val());
  EXPECT_FLOAT_EQ(
      0.54899865,
      stan::math::multi_normal_cholesky_lpdf(y, mu, L).d_.val_.val());
}
