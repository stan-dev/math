#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/fwd/scal.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, binormal_integral_using) {
  using stan::math::std_binormal_integral;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan_fd) {
  using stan::math::fvar;
  fvar<double> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<double> rho = 0.3;
  fvar<double> a = nan;
  fvar<double> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_fd) {
  using stan::math::fvar;
  fvar<double> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<double> rho = 0.3;
  fvar<double> a = 2;
  fvar<double> b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_fd) {
  using stan::math::fvar;
  fvar<double> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<double> rho = nan;
  fvar<double> a = 2;
  fvar<double> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_fd) {
  using stan::math::fvar;
  fvar<double> rho = -1.3;
  fvar<double> a = 2;
  fvar<double> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_fd) {
  using stan::math::fvar;
  fvar<double> rho = 1.3;
  fvar<double> a = 2;
  fvar<double> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_fd) {
  using stan::math::fvar;
  fvar<double> rho = 0.3;
  fvar<double> a = 2;
  fvar<double> b = 1;
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
//// ffv
TEST(MathFunctions, binormal_integral_throw_RV_1_nan_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = nan;
  fvar<fvar<double>> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<double>> rho = nan;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = -1.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = 1.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_test_fd) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  using stan::math::fvar;
  fvar<double> rho = 0.3;
  fvar<double> a = -0.4;
  fvar<double> b = 2.7;
  fvar<double> f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.344276561500873, f.val_);

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.3445782583896758, f.val_);

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9937227710497979, f.val_);

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9997643606337163, f.val_);

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.146032113348184e-06, f.val_);

  rho = -0.99;
  a = -4.5;
  b = 10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(3.397673124738709e-06, f.val_);

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(4.016000583859118e-11, f.val_);

  rho = -0.99;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_NEAR(0, f.val_, 1e-16);

  rho = 0.99;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(7.619853024160583e-24, f.val_);

  rho = 0.5;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(5.612932952882069e-24, f.val_);

  rho = 0.5;
  a = -4.5;
  b = -4.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.348555308541017e-08, f.val_);

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val_);

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val_);

  rho = -0.5;
  a = -3.5;
  b = -3.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(8.061334595291052e-14, f.val_);

  rho = -0.3;
  a = 2.3;
  b = -3.4;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(std::log(0.0003017720314683548), std::log(f.val_));
}

TEST(MathFunctions, binormal_integral_val_test_ffd) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  using stan::math::fvar;
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = -0.4;
  fvar<fvar<double>> b = 2.7;
  fvar<fvar<double>> f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.344276561500873, f.val_.val_);

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.3445782583896758, f.val_.val_);

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9937227710497979, f.val_.val_);

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(0.9997643606337163, f.val_.val_);

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.146032113348184e-06, f.val_.val_);

  rho = -0.99;
  a = -4.5;
  b = 10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(3.397673124738709e-06, f.val_.val_);

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(4.016000583859118e-11, f.val_.val_);

  rho = -0.99;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_NEAR(0, f.val_.val_, 1e-16);

  rho = 0.99;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(7.619853024160583e-24, f.val_.val_);

  rho = 0.5;
  a = -4.5;
  b = -10;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(5.612932952882069e-24, f.val_.val_);

  rho = 0.5;
  a = -4.5;
  b = -4.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.348555308541017e-08, f.val_.val_);

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val_.val_);

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val_.val_);

  rho = -0.5;
  a = -3.5;
  b = -3.5;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(8.061334595291052e-14, f.val_.val_);

  rho = -0.3;
  a = 2.3;
  b = -3.4;
  f = stan::math::std_binormal_integral(a, b, rho);
  EXPECT_FLOAT_EQ(std::log(0.0003017720314683548), std::log(f.val_.val_));
}
