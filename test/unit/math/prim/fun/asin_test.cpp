#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, asin) {
  stan::test::expect_common_prim([](auto x) { return std::asin(x); },
                                 [](auto x) { return stan::math::asin(x); });
}
