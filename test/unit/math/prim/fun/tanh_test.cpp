#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, tanh) {
  stan::test::expect_common_prim([](auto x) { return std::tanh(x); },
                                 [](auto x) { return stan::math::tanh(x); });
}
