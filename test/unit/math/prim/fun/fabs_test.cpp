#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, fabs) {
  stan::test::expect_common_prim([](auto x) { return std::fabs(x); },
                                 [](auto x) { return stan::math::fabs(x); });
}
