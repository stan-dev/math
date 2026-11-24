#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, logInvLogitDiff) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_inv_logit_diff(x1, x2);
  };
  stan::test::expect_ad(f, 0.5, -1.0);
  stan::test::expect_ad(f, 0.5, 0.0);
  stan::test::expect_ad(f, 1.2, 0.6);
  stan::test::expect_ad(f, 3.4, 0.9);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}

TEST(mathMixScalFun, logInvLogitDiff_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_inv_logit_diff;
    return log_inv_logit_diff(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, -3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

TEST(mathMixScalFun, logInvLogitDiffGradientLargeArgs) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;
  // {x, y, d/dx, d/dy}: mpmath, 50 digits
  std::vector<std::vector<double>> points{
      {40, 0, 8.4967085105831779907e-18, -5.0000000000000000425e-1},
      {100, 1, 1.3832290902125321262e-43, -7.3105857863000487925e-1},
      {30, -30, 9.3576229688401746049e-14, -9.3576229688401746049e-14},
      {37, 36, 5.8197670686932650972e-1, -1.5819767068693261924},
      {-30, -40, 1.0000454019909161115, -4.5401991009692016683e-5},
      {1, 0, 8.5091812823932154513e-1, -1.0819767068693264244},
      {700, 20, 4.7835719068902113962e-296, -9.9999999793884638181e-1},
      {0.5, -1.0, 6.647575855870136797e-1, -5.5615833815886336509e-1}};
  for (const auto& p : points) {
    double x = p[0];
    double y = p[1];
    double dx = p[2];
    double dy = p[3];

    var xv = x;
    var yv = y;
    log_inv_logit_diff(xv, yv).grad();
    EXPECT_NEAR(xv.adj(), dx, 1e-13 * dx) << x << ", " << y;
    EXPECT_NEAR(yv.adj(), dy, -1e-13 * dy) << x << ", " << y;
    stan::math::recover_memory();

    var xv2 = x;
    log_inv_logit_diff(xv2, y).grad();
    EXPECT_NEAR(xv2.adj(), dx, 1e-13 * dx) << x << ", " << y;
    stan::math::recover_memory();

    EXPECT_NEAR(log_inv_logit_diff(fvar<double>(x, 1), y).d_, dx, 1e-13 * dx)
        << x << ", " << y;
    EXPECT_NEAR(log_inv_logit_diff(fvar<double>(x, 1), fvar<double>(y, 0)).d_,
                dx, 1e-13 * dx)
        << x << ", " << y;
    EXPECT_NEAR(log_inv_logit_diff(fvar<double>(x, 0), fvar<double>(y, 1)).d_,
                dy, -1e-13 * dy)
        << x << ", " << y;
  }
}

TEST(mathMixScalFun, logInvLogitDiffGradientEqualArgs) {
  using stan::math::fvar;
  using stan::math::log_inv_logit_diff;
  using stan::math::var;
  // behaviour change: d/dx at x == y is +inf, the limit from x > y (was -inf)
  double inf = std::numeric_limits<double>::infinity();
  var xv = 2.0;
  var yv = 2.0;
  log_inv_logit_diff(xv, yv).grad();
  EXPECT_EQ(xv.adj(), inf);
  EXPECT_EQ(yv.adj(), -inf);
  stan::math::recover_memory();

  var xv2 = 2.0;
  log_inv_logit_diff(xv2, 2.0).grad();
  EXPECT_EQ(xv2.adj(), inf);
  stan::math::recover_memory();

  EXPECT_EQ(log_inv_logit_diff(fvar<double>(2, 1), 2.0).d_, inf);
}
