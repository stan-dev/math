#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, imag) {
  auto f = [](const auto& z) { return imag(z); };
  stan::test::expect_complex_common(f);
}
