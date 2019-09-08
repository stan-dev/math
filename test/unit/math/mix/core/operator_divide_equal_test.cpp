#include <test/unit/math/test_ad.hpp>

TEST(mathMixCore, operatorDivideEqual) {
  auto f = [](const auto& x1, const auto& x2) {
    // decltype instead of auto to make following statement legal
    decltype(x1 + x2) y = x1;
    y /= x2;
    return y;
  };
  bool is_assignment = true;
  stan::test::expect_common_binary(f, is_assignment);
}
