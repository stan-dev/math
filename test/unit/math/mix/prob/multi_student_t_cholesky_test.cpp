#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsMultiStudentTCholesky, matvar) {
  auto f = [](const auto& y, const auto& nu, const auto& mu, const auto& L) {
    return stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L);
  };

  auto f_const_y = [](const auto& y) {
    return [&y](const auto& nu, const auto& mu, const auto& L) {
      return stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L);
    };
  };

  auto f_const_nu = [](const auto& nu) {
    return [&nu](const auto& y, const auto& mu, const auto& L) {
      return stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L);
    };
  };

  Eigen::VectorXd y1(1);
  y1 << 1;
  Eigen::VectorXd mu1(1);
  mu1 << 3.4;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  double nu = 1.1;

  stan::test::expect_ad(f_const_y(y1), nu, mu1, Sigma11);
  stan::test::expect_ad(f_const_nu(nu), y1, mu1, Sigma11);
  stan::test::expect_ad_matvar(f, y1, nu, mu1, Sigma11);

  Eigen::VectorXd y0(0);
  Eigen::VectorXd mu0(0);
  Eigen::MatrixXd Sigma00(0, 0);
  stan::test::expect_ad(f_const_y(y0), nu, mu0, Sigma00);
  stan::test::expect_ad(f_const_nu(nu), y0, mu0, Sigma00);
  stan::test::expect_ad_matvar(f, y0, nu, mu0, Sigma00);

  Eigen::VectorXd y2(2);
  y2 << 1.0, 0.1;
  Eigen::VectorXd mu2(2);
  mu2 << 0.1, 2.0;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 2.0, 0.5, 0.5, 1.1;
  stan::test::expect_ad(f_const_y(y2), nu, mu2, Sigma22);
  stan::test::expect_ad(f_const_nu(nu), y2, mu2, Sigma22);
  stan::test::expect_ad_matvar(f, y2, nu, mu2, Sigma22);

  Eigen::VectorXd y22(2);
  y22 << 0.4, 0.3;
  Eigen::VectorXd mu22(2);
  mu22 << 2.1, 1.0;
  std::vector<Eigen::VectorXd> y2s = {y2, y22};
  std::vector<Eigen::VectorXd> mu2s = {mu2, mu22};
  stan::test::expect_ad(f_const_y(y2), nu, mu2, Sigma22);
  stan::test::expect_ad(f_const_y(y2), nu, mu2s, Sigma22);
  stan::test::expect_ad(f_const_y(y2s), nu, mu2s, Sigma22);
  stan::test::expect_ad(f_const_nu(nu), y2, mu2, Sigma22);
  stan::test::expect_ad(f_const_nu(nu), y2s, mu2, Sigma22);
  stan::test::expect_ad(f_const_nu(nu), y2s, mu2s, Sigma22);
  stan::test::expect_ad(f_const_nu(nu), y2, mu2s, Sigma22);
  stan::test::expect_ad_matvar(f, y2s, nu, mu2s, Sigma22);
  stan::test::expect_ad_matvar(f, y2, nu, mu2s, Sigma22);
  stan::test::expect_ad_matvar(f, y2s, nu, mu2, Sigma22);

  // Error sizes
  stan::test::expect_ad(f_const_y(y0), nu, mu0, Sigma11);
  stan::test::expect_ad(f_const_nu(nu), y0, mu0, Sigma11);
  stan::test::expect_ad(f_const_y(y1), nu, mu0, Sigma11);
  stan::test::expect_ad(f_const_nu(nu), y1, mu0, Sigma11);
  stan::test::expect_ad_matvar(f, y0, nu, mu0, Sigma11);
  stan::test::expect_ad_matvar(f, y1, nu, mu1, Sigma00);
}

TEST(ProbDistributionsMultiStudentTCholesky, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::multi_student_t_log;
  using stan::math::var;
  using std::vector;
  Matrix<fvar<var>, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<fvar<var>, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  double nu = 4.0;

  for (int i = 0; i < 3; i++) {
    y(i).d_ = 1.0;
    mu(i).d_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_ = 1.0;
  }

  Matrix<fvar<var>, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  fvar<var> lp = multi_student_t_cholesky_lpdf(y, nu, mu, L);
  EXPECT_NEAR(-10.1246, lp.val_.val(), 0.0001);
  EXPECT_NEAR(-0.0411685, lp.d_.val(), 0.0001);

  stan::math::recover_memory();
}

TEST(ProbDistributionsMultiStudentTCholesky, fvar_fvar_var) {
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
  double nu = 4.0;

  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  for (int i = 0; i < 3; i++) {
    y(i).d_.val_ = 1.0;
    mu(i).d_.val_ = 1.0;
    for (int j = 0; j < 3; j++)
      Sigma(i, j).d_.val_ = 1.0;
  }

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  EXPECT_NEAR(
      -10.1246,
      stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L).val_.val_.val(),
      0.0001);
  EXPECT_NEAR(
      -0.0411685,
      stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L).d_.val_.val(),
      0.0001);
}
