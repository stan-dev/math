#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix,  iTimes) {
  auto f = [](const auto& x) { return stan::math::i_times(x); };
  stan::test::expect_complex_common(f);
}

TEST_F(mathMix,  negITimes) {
  auto f = [](const auto& x) { return stan::math::neg_i_times(x); };
  stan::test::expect_complex_common(f);
}

TEST_F(mathMix,  complexNegate) {
  auto f = [](const auto& x) { return stan::math::complex_negate(x); };
  stan::test::expect_complex_common(f);
}
