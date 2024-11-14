#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, operatorGreaterThanOrEqual) {
  auto f = [](const auto& x1, const auto& x2) { return x1 >= x2; };
  stan::test::expect_common_comparison(f);
}
