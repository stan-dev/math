#include <test/unit/math/test_ad.hpp>
#include <array>
#include <cmath>
#include <vector>

TEST(mathMixMatFun, tanh) {
  auto f = [](const auto& x1) {
    using stan::math::tanh;
    return tanh(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized<
      stan::test::ScalarSupport::RealAndComplex>(f, -2.6, -2, -1.2, -0.5, 0.5,
                                                 1.5);
  stan::test::expect_complex_common(f);
}

// First and second derivatives at large |x| against mpmath. 1 - tanh(x)^2
// is 0 in double from |x| = 19.06, and 1 / cosh(x)^2 makes the second
// derivative NaN once cosh(x)^2 overflows (|x| > 355).
TEST(mathMixMatFun, tanh_derivatives_large_argument) {
  using stan::math::fvar;
  using stan::math::var;
  // x, d/dx tanh(x), d2/dx2 tanh(x)
  std::vector<std::array<double, 3>> cases
      = {{3.0, 0.0098660371654401913, -0.019634494363042436},
         {-20.0, 1.6993417021166356e-17, 3.3986834042332711e-17},
         {20.0, 1.6993417021166356e-17, -3.3986834042332711e-17},
         {300.0, 1.0601586212017243e-260, -2.1203172424034487e-260},
         {1000.0, 0.0, 0.0}};
  for (const auto& c : cases) {
    const double x = c[0];
    const double d1 = c[1];
    const double d2 = c[2];

    fvar<double> y1 = stan::math::tanh(fvar<double>(x, 1.0));
    EXPECT_NEAR(d1, y1.d_, 1e-14 * std::fabs(d1)) << "x = " << x;

    fvar<fvar<double>> y2 = stan::math::tanh(
        fvar<fvar<double>>(fvar<double>(x, 1.0), fvar<double>(1.0, 0.0)));
    EXPECT_NEAR(d1, y2.d_.val_, 1e-14 * std::fabs(d1)) << "x = " << x;
    EXPECT_NEAR(d2, y2.d_.d_, 1e-14 * std::fabs(d2)) << "x = " << x;

    fvar<var> x3(var(x), 1.0);
    fvar<var> y3 = stan::math::tanh(x3);
    y3.d_.grad();
    EXPECT_NEAR(d1, y3.d_.val(), 1e-14 * std::fabs(d1)) << "x = " << x;
    EXPECT_NEAR(d2, x3.val_.adj(), 1e-14 * std::fabs(d2)) << "x = " << x;
    stan::math::recover_memory();
  }
}

// First, second and third derivatives in fvar<fvar<var>> against mpmath,
// at 0 (where an fabs(x) in the derivative would make the third derivative
// 0), on both sides of 0.5, and in the tail.
TEST(mathMixMatFun, tanh_third_derivative) {
  using stan::math::fvar;
  using stan::math::var;
  // x, tanh'(x), tanh''(x), tanh'''(x)
  std::vector<std::array<double, 4>> cases = {
      {0.0, 1.0, 0.0, -2.0},
      {0.25, 0.94001484880637796, -0.46045435881856592, -1.5417081006333538},
      {0.4999, 0.78652041633736061, -0.72680544069311918, -0.56560452654352799},
      {0.5, 0.78644773296592741, -0.72686198138358728, -0.56520928825977036},
      {-0.5, 0.78644773296592741, 0.72686198138358728, -0.56520928825977036},
      {3.0, 0.0098660371654401913, -0.019634494363042436, 0.038880116525661682},
      {-20.0, 1.6993417021166356e-17, 3.3986834042332711e-17,
       6.7973668084665422e-17}};
  for (const auto& c : cases) {
    const double x = c[0];
    var x_v = x;
    fvar<fvar<var>> x_ffv(fvar<var>(x_v, 1.0), fvar<var>(1.0, 0.0));
    fvar<fvar<var>> y = stan::math::tanh(x_ffv);
    y.d_.d_.grad();
    EXPECT_NEAR(c[1], y.d_.val_.val(), 1e-14 * std::fabs(c[1])) << "x = " << x;
    EXPECT_NEAR(c[2], y.d_.d_.val(), 1e-14 * std::fabs(c[2]) + 1e-300)
        << "x = " << x;
    EXPECT_NEAR(c[3], x_v.adj(), 1e-14 * std::fabs(c[3])) << "x = " << x;
    stan::math::recover_memory();
  }
}

TEST(mathMixMatFun, tanh_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::tanh;
    return tanh(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-2.6, -0.5, 0.5, 1.5};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
