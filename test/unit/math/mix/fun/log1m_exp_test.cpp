#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1m_exp) {
  auto f = [](const auto& x1) { return stan::math::log1m_exp(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -14, -12.6, -2, -1, -0.2, -0.5, 1.3,
                                      3);
}
