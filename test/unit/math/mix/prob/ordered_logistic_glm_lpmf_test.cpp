#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_ordered_logistic_glm_lpmf) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& beta, const auto& cutpoints) {
      return stan::math::ordered_logistic_glm_lpmf(y, x, beta, cutpoints);
    };
  };

  std::vector<int> y{0, 1, 2};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd cutpoints = stan::math::sort_asc(Eigen::VectorXd::Random(2));
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y[0]), x, beta, cutpoints);
  stan::test::expect_ad(f(y[0]), x_rowvec, beta, cutpoints);
  stan::test::expect_ad(f(y), x, beta, cutpoints);
  stan::test::expect_ad(f(y), x_rowvec, beta, cutpoints);
}

// Regression test for the higher-order-AD NaN bug reported 2026-05-21
// (mirror of the test in ordered_logistic_test.cpp; same root cause:
// orphan exp(+INF) vari created by the partials block for y == N_classes
// and y == 1 sentinel cases). See that file for the full bug description.
TEST_F(AgradRev,
       mathMixScalFun_ordered_logistic_glm_lpmf_top_category_higher_order_ad) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_glm_lpmf;
  using stan::math::var;

  // N_classes = 3 (cutpoints size 2). y contains both 1 and N_classes.
  std::vector<int> y{1, 2, 3, 2, 1};
  const int N = 5;
  const int D = 2;

  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, Eigen::Dynamic> x_ffv(N, D);
  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> beta_ffv(D);
  Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> cuts_ffv(2);

  x_ffv << 0.5, -0.3, 0.1, 0.7, -0.2, 0.4, 0.8, -0.1, 0.3, 0.2;
  beta_ffv << 0.4, -0.6;
  cuts_ffv << 0.0, 1.0;

  // Mixed-zero tangent seeds (compute_s2 / third_diff pattern).
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < D; j++) {
      x_ffv(i, j).d_ = (i + j) % 2 == 0 ? 1.0 : 0.0;
      x_ffv(i, j).val_.d_ = (i + j) % 3 == 0 ? 1.0 : 0.0;
    }
  }
  for (int j = 0; j < D; j++) {
    beta_ffv(j).d_ = 1.0;
    beta_ffv(j).val_.d_ = 0.0;
  }
  for (int j = 0; j < 2; j++) {
    cuts_ffv(j).d_ = 0.0;
    cuts_ffv(j).val_.d_ = 1.0;
  }

  fvar<fvar<var>> out_ffv
      = ordered_logistic_glm_lpmf(y, x_ffv, beta_ffv, cuts_ffv);

  EXPECT_TRUE(std::isfinite(out_ffv.val_.val_.val()));
  EXPECT_TRUE(std::isfinite(out_ffv.d_.val_.val()));
  EXPECT_TRUE(std::isfinite(out_ffv.val_.d_.val()));
  EXPECT_TRUE(std::isfinite(out_ffv.d_.d_.val()));

  out_ffv.d_.d_.grad();

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < D; j++) {
      EXPECT_TRUE(std::isfinite(x_ffv(i, j).val_.val_.adj()))
          << "x_ffv(" << i << "," << j
          << ").val_.val_.adj() non-finite (y[" << i << "]=" << y[i] << ")";
    }
  }
  for (int j = 0; j < D; j++) {
    EXPECT_TRUE(std::isfinite(beta_ffv(j).val_.val_.adj()))
        << "beta_ffv(" << j << ").val_.val_.adj() non-finite";
  }
  for (int j = 0; j < 2; j++) {
    EXPECT_TRUE(std::isfinite(cuts_ffv(j).val_.val_.adj()))
        << "cuts_ffv(" << j << ").val_.val_.adj() non-finite";
  }
}
