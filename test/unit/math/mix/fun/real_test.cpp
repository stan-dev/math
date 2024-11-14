#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, real) {
  auto f = [](const auto& z) { return real(z); };
  stan::test::expect_complex_common(f);
}
