#include <test/unit/math/test_ad.hpp>

TEST(matPrimFun, iTimes) {
  auto f = [](const auto& x) { return stan::math::i_times(x); };
  stan::test::expect_complex_common(f);
}

TEST(matPrimFun, negITimes) {
  auto f = [](const auto& x) { return stan::math::neg_i_times(x); };
  stan::test::expect_complex_common(f);
}

TEST(matPrimFun, complexNegate) {
  auto f = [](const auto& x) { return stan::math::complex_negate(x); };
  stan::test::expect_complex_common(f);
}
