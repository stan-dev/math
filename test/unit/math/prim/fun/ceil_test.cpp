#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, ceil) {
  stan::test::expect_common_prim([](auto x) { return std::ceil(x); },
                                 [](auto x) { return stan::math::ceil(x); });
}
