#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, opratorPlusPlusPre1) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = ++y;
    return y;
  };
  stan::test::expect_common_unary(f);
}
TEST(mathMixCore, opratorPlusPlusPre2) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = ++y;
    return z;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, opratorPlusPlusPost1) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y++;
    return y;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, opratorPlusPlusPost2) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y++;
    return z;
  };
  stan::test::expect_common_unary(f);
}
