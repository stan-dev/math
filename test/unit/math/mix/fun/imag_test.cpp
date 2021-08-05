#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, imag) {
  auto f = [](const auto& z) { return stan::math::imag(z); };
  stan::test::expect_complex_common(f);
}
