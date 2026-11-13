#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, besselI0eDerivativeRatio) {
  auto f = [](const auto& x) {
    return stan::math::bessel_i0e_derivative_ratio(0, x);
  };
  stan::test::expect_unary_vectorized(f, 0.5, 2.0, 7.7, 15.0, 30.0, 100.0);
  stan::test::expect_ad(f, -1.0);
  auto f2 = [](const auto& x) {
    return stan::math::bessel_i0e_derivative_ratio(2, x);
  };
  stan::test::expect_ad(f2, 0.5);
  stan::test::expect_ad(f2, 15.0);
  // the largest order is accepted for autodiff arguments too
  stan::math::var z = 480.0;
  EXPECT_NO_THROW(stan::math::bessel_i0e_derivative_ratio(1000, z));
  stan::math::recover_memory();
}

TEST(mathMixScalFun, besselI0eDerivativeRatioThirdOrder) {
  using stan::math::fvar;
  using stan::math::var;
  // z, then B = 1 - I1(z) / I0(z) and its first three derivatives, computed
  // with mpmath (60 digits); finite differences cannot check large z
  const double ref[][5] = {{0.5, 0.75750038741919805, -0.45619471273655707,
                            0.16364505734995701, 0.24000549054178479},
                           {15, 0.033930436013491875, -0.0023049599398175424,
                            0.00031352477751293574, -6.4052326658463717e-5},
                           {1e6, 5.00000125000125e-7, -5.00000250000375e-13,
                            1.0000007500015e-18, -3.0000030000075e-24},
                           {1e12, 5.00000000000125e-13, -5.0000000000025e-25,
                            1.00000000000075e-36, -3.000000000003e-48}};
  for (const auto& r : ref) {
    fvar<fvar<var>> z(fvar<var>(r[0], 1.0), fvar<var>(1.0, 0.0));
    fvar<fvar<var>> b = stan::math::bessel_i0e_derivative_ratio(0, z);
    // below 25, the value comes from approximations of I0 and I1
    const double tol = r[0] < 25 ? 1e-13 : 1e-14;
    EXPECT_NEAR(r[1], b.val_.val_.val(), tol * std::fabs(r[1]));
    EXPECT_NEAR(r[2], b.d_.val_.val(), 1e-13 * std::fabs(r[2]));
    EXPECT_NEAR(r[3], b.d_.d_.val(), 1e-13 * std::fabs(r[3]));
    b.d_.d_.grad();
    EXPECT_NEAR(r[4], z.val_.val_.adj(), 1e-12 * std::fabs(r[4]));
    stan::math::recover_memory();
  }
}
