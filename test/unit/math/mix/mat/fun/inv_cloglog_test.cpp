#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invCLogLog) {
  auto f = [](const auto& x1) { return stan::math::inv_cloglog(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
