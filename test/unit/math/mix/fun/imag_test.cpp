#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, imag) {
  auto f = [](const auto& z) { return imag(z); };
  stan::test::expect_complex_common(f);
}
