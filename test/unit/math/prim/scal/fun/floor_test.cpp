#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, floor) {
  stan::test::expect_common_prim([](auto x) { return std::floor(x); },
                                 [](auto x) { return stan::math::floor(x); });
}
