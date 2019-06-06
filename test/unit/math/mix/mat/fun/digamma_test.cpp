#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, digamma) {
  auto f = [](const auto& x1) { return stan::math::digamma(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
}
