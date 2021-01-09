#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, step) {
  stan::test::expect_common_prim([](auto x) { return x < 0.0 ? 0 : 1; },
                                 [](auto x) { return stan::math::step(x); });
}
