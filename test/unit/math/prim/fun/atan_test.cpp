#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, atan) {
  stan::test::expect_common_prim([](auto x) { return std::atan(x); },
                                 [](auto x) { return stan::math::atan(x); });
}
