#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, real) {
  auto f = [](const auto& z) { return stan::math::real(z); };
  stan::test::expect_complex_common(f);
}
