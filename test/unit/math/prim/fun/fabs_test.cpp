#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, fabs) {
  stan::test::expect_common_prim([](auto x) { return std::fabs(x); },
                                 [](auto x) { return stan::math::fabs(x); });
}

TEST(primScalFun, fabs_complex) {
  using stan::math::fabs;
  std::complex<double> z(1.0, 2.0);
  EXPECT_FLOAT_EQ(fabs(z), 2.236068);
}
