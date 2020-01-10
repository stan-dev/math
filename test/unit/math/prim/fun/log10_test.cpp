#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, log10) {
  stan::test::expect_common_prim([](auto x) { return std::log10(x); },
                                 [](auto x) { return stan::math::log10(x); });
}
