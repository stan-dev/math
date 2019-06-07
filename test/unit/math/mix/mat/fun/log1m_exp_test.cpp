#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1m_exp) {
  auto f = [](const auto& x1) { return stan::math::log1m_exp(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
