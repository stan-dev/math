#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, sinh) {
  stan::test::expect_common_prim([](auto x) { return std::sinh(x); },
                                 [](auto x) { return stan::math::sinh(x); });
}
