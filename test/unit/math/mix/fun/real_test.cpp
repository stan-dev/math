#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, real) {
  auto f = [](const auto& z) { return real(z); };
  stan::test::expect_complex_common(f);
}
