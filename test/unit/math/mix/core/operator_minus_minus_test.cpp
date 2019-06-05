#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, opratorMinusMinusPre1) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = --y;
    return y;
  };
  stan::test::expect_common_unary(f);
}
TEST(mathMixCore, opratorMinusMinusPre2) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = --y;
    return z;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, opratorMinusMinusPost1) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y--;
    return y;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, opratorMinusMinusPost2) {
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y--;
    return z;
  };
  stan::test::expect_common_unary(f);
}
