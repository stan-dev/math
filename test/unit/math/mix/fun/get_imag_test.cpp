#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, get_imag) {
  auto f = [](const auto& z) { return stan::math::get_imag(z); };
  stan::test::expect_complex_common(f);
}
