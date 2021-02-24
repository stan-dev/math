#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsInvWishart, matvar) {
  auto f = [](const auto& y, const auto& dof, const auto& sigma) {
    auto y_sym = stan::math::multiply(0.5, y + y.transpose());
    auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
    return stan::math::inv_wishart_lpdf(y_sym, dof, sigma_sym);
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

TEST(ProbDistributionsInvWishart, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::inv_wishart_log;
  using stan::math::var;

  Matrix<fvar<var>, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.0910458, -11.9036079, 16.7585782, 0.8530256,
      1.0910458, 0.8530256, 2.5786609;

  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.0597878, 1.107166, 3.0597878, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i, j).d_ = 1.0;
      Sigma(i, j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_log(Y, dof, Sigma).val_.val(),
              0.01);
  EXPECT_NEAR(-1.4893348387330674,
              stan::math::inv_wishart_log(Y, dof, Sigma).d_.val(), 0.01);

  stan::math::recover_memory();
}

TEST(ProbDistributionsInvWishart, fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::inv_wishart_log;
  using stan::math::var;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Y(3, 3);
  Y << 12.147233, -11.9036079, 1.0910458, -11.9036079, 16.7585782, 0.8530256,
      1.0910458, 0.8530256, 2.5786609;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 7.785215, 3.059788, 1.107166, 3.059788, 10.3515035, -0.1232598,
      1.107166, -0.1232598, 7.7623386;

  double dof = 4.0;
  double log_p = log(2.008407e-08);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      Y(i, j).d_ = 1.0;
      Sigma(i, j).d_ = 1.0;
    }

  EXPECT_NEAR(log_p, stan::math::inv_wishart_log(Y, dof, Sigma).val_.val_.val(),
              0.01);
  EXPECT_NEAR(-1.4893348387330674,
              stan::math::inv_wishart_log(Y, dof, Sigma).d_.val_.val(), 0.01);

  stan::math::recover_memory();
}
