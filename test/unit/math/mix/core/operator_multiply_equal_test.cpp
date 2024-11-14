#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, operatorMultiplyEqual) {
  auto f = [](const auto& x1, const auto& x2) {
    decltype(x1 + x2) y = x1;
    y *= x2;
    return y;
  };
  bool is_assignment = true;
  stan::test::expect_common_binary(f, is_assignment);
}
