#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, cos) {
  stan::test::expect_common_prim([](auto x) { return std::cos(x); },
                                 [](auto x) { return stan::math::cos(x); });
}
