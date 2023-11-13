#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorPlusPlusPre1) {
  // this functor tests that y is right value after operator
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = ++y;
    // Suppress unused variable warning
    z = z + 0;
    return y;
  };
  stan::test::expect_common_unary(f);
}
TEST(mathMixCore, operatorPlusPlusPre2) {
  // this functor tests that value of expression has right value
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = ++y;
    return z;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, operatorPlusPlusPost1) {
  // this functor tests that y is right value after operator
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y++;
    // Suppress unused variable warning
    z = z + 0;
    return y;
  };
  stan::test::expect_common_unary(f);
}

TEST(mathMixCore, operatorPlusPlusPost2) {
  // this functor tests that value of expression has right value
  auto f = [](const auto& x1) {
    auto y = x1;
    auto z = y++;
    return z;
  };
  stan::test::expect_common_unary(f);
}
