#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, to_complex) {
  auto f = [](const auto& x) { return stan::math::to_complex(x); };

  stan::test::expect_common_unary(f);
}
