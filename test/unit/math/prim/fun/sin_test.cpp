#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, sin) {
  stan::test::expect_common_prim([](auto x) { return std::sin(x); },
                                 [](auto x) { return stan::math::sin(x); });
}
