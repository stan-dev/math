#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, logBesselI0e) {
  auto f = [](const auto& x) { return stan::math::log_bessel_i0e(x); };
  stan::test::expect_unary_vectorized(f, 0.5, 2.0, 7.7, 7.8, 30.0, 400.0,
                                      600.0);
  stan::test::expect_ad(f, -1.0);
}

TEST(mathMixScalFun, logBesselI0eThirdOrder) {
  using stan::math::fvar;
  using stan::math::var;
  // z, B = 1 - I1(z) / I0(z), B' and B'' from mpmath (60 digits); the first
  // three derivatives of log(I0(z)) - z are -B, -B' and -B''
  const double ref[][4] = {
      {0.5, 0.75750038741919805, -0.45619471273655707, 0.16364505734995701},
      {15, 0.033930436013491875, -0.0023049599398175424,
       0.00031352477751293574},
      {1e6, 5.00000125000125e-7, -5.00000250000375e-13, 1.0000007500015e-18},
      {1e12, 5.00000000000125e-13, -5.0000000000025e-25, 1.00000000000075e-36}};
  for (const auto& r : ref) {
    fvar<fvar<var>> z(fvar<var>(r[0], 1.0), fvar<var>(1.0, 0.0));
    fvar<fvar<var>> l = stan::math::log_bessel_i0e(z);
    // below 25, B comes from approximations of I0 and I1
    const double tol = r[0] < 25 ? 1e-13 : 1e-14;
    EXPECT_NEAR(-r[1], l.d_.val_.val(), tol * std::fabs(r[1]));
    EXPECT_NEAR(-r[1], l.val_.d_.val(), tol * std::fabs(r[1]));
    EXPECT_NEAR(-r[2], l.d_.d_.val(), 1e-13 * std::fabs(r[2]));
    l.d_.d_.grad();
    EXPECT_NEAR(-r[3], z.val_.val_.adj(), 1e-13 * std::fabs(r[3]));
    stan::math::recover_memory();
  }
}
