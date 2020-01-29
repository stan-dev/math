#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, acos) {
  stan::test::expect_common_prim([](auto x) { return std::acos(x); },
                                 [](auto x) { return stan::math::acos(x); });
}
