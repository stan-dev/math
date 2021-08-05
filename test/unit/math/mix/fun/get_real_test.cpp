#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, get_real) {
  auto f = [](const auto& z) { return stan::math::get_real(z); };
  stan::test::expect_complex_common(f);
}
