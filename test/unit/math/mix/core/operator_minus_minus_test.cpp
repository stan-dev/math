#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorMinusMinusPre1) {
  // this functor tests that y is right value after operator
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = --y;
    // Suppress unused variable warning
    z = z + 0;
    return y;
  };
  stan::test::expect_common_unary(f);
}
TEST(mathMixCore, operatorMinusMinusPre2) {
  // this functor tests that value of expression has right value
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = --y;
    return z;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, operatorMinusMinusPost1) {
  // this functor tests that y is right value after operator
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y--;
    // Suppress unused variable warning
    z = z + 0;
    return y;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, operatorMinusMinusPost2) {
  // this functor tests that value of expression has right value
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y--;
    return z;
  };
  stan::test::expect_common_unary(f);
}
