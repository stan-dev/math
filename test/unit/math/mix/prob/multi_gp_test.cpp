#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsMultiGP, matvar) {
  auto f = [](const auto& y, const auto& sigma, const auto& w) {
    auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
    return stan::math::multi_gp_lpdf(y, sigma_sym, w);
  };

  Eigen::MatrixXd y11(1, 1);
  y11 << 1;
  Eigen::VectorXd w1(1);
  w1 << 3.4;
  Eigen::MatrixXd Sigma11(1, 1);
  Sigma11 << 1;
  stan::test::expect_ad(f, y11, Sigma11, w1);
  stan::test::expect_ad_matvar(f, y11, Sigma11, w1);

  Eigen::MatrixXd y00(0, 0);
  Eigen::VectorXd w0(0);
  Eigen::MatrixXd Sigma00(0, 0);
  stan::test::expect_ad(f, y00, Sigma00, w0);
  stan::test::expect_ad_matvar(f, y00, Sigma00, w0);

  Eigen::MatrixXd y22(2, 2);
  y22 << 1.0, 0.1, 0.3, 0.5;
  Eigen::VectorXd w2(2);
  w2 << 0.1, 2.0;
  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 3.0, 0.5, 0.5, 2.0;
  stan::test::expect_ad(f, y22, Sigma22, w2);
  stan::test::expect_ad_matvar(f, y22, Sigma22, w2);

  Eigen::MatrixXd y12(1, 2);
  y12 << 1.0, 0.1;
  stan::test::expect_ad(f, y12, Sigma22, w1);
  stan::test::expect_ad_matvar(f, y12, Sigma22, w1);

  // Error sizes
  stan::test::expect_ad(f, y00, Sigma00, w2);
  stan::test::expect_ad(f, y00, Sigma22, w0);
  stan::test::expect_ad(f, y00, Sigma22, w2);
  stan::test::expect_ad(f, y22, Sigma22, w0);
  stan::test::expect_ad(f, y22, Sigma00, w2);
  stan::test::expect_ad(f, y22, Sigma00, w0);
}

TEST(ProbDistributionsMultiGP, fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<var>, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<var>, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<var>, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<fvar<var>, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  for (int i = 0; i < 5; i++) {
    mu(i).d_ = 1.0;
    if (i < 3)
      w(i).d_ = 1.0;
    for (int j = 0; j < 5; j++) {
      Sigma(i, j).d_ = 1.0;
      if (i < 3)
        y(i, j).d_ = 1.0;
    }
  }

  fvar<var> lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<var>, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<var>, Dynamic, Dynamic> cSigma((1.0 / w[i]) * Sigma);
    lp_ref += stan::math::multi_normal_log(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_.val(),
                  stan::math::multi_gp_log(y, Sigma, w).val_.val());
  EXPECT_FLOAT_EQ(-74.572952, stan::math::multi_gp_log(y, Sigma, w).d_.val());

  stan::math::recover_memory();
}

TEST(ProbDistributionsMultiGP, fvar_fvar_var) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::var;
  Matrix<fvar<fvar<var> >, Dynamic, 1> mu(5, 1);
  mu.setZero();

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> y(3, 5);
  y << 2.0, -2.0, 11.0, 4.0, -2.0, 11.0, 2.0, -5.0, 11.0, 0.0, -2.0, 11.0, 2.0,
      -2.0, -11.0;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> Sigma(5, 5);
  Sigma << 9.0, -3.0, 0.0, 0.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0,
      1.0, 0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0;

  Matrix<fvar<fvar<var> >, Dynamic, 1> w(3, 1);
  w << 1.0, 0.5, 1.5;

  for (int i = 0; i < 5; i++) {
    mu(i).d_.val_ = 1.0;
    if (i < 3)
      w(i).d_.val_ = 1.0;
    for (int j = 0; j < 5; j++) {
      Sigma(i, j).d_.val_ = 1.0;
      if (i < 3)
        y(i, j).d_.val_ = 1.0;
    }
  }

  fvar<fvar<var> > lp_ref(0);
  for (size_t i = 0; i < 3; i++) {
    Matrix<fvar<fvar<var> >, Dynamic, 1> cy(y.row(i).transpose());
    Matrix<fvar<fvar<var> >, Dynamic, Dynamic> cSigma((1.0 / w[i]) * Sigma);
    lp_ref += stan::math::multi_normal_log(cy, mu, cSigma);
  }

  EXPECT_FLOAT_EQ(lp_ref.val_.val_.val(),
                  stan::math::multi_gp_log(y, Sigma, w).val_.val_.val());
  EXPECT_FLOAT_EQ(-74.572952,
                  stan::math::multi_gp_log(y, Sigma, w).d_.val_.val());

  stan::math::recover_memory();
}
