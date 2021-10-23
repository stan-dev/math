#include <test/unit/math/test_ad.hpp>
#include <stan/math/rev/fun/lmultiply.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, lmultiply) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::lmultiply(x1, x2);
  };
  for (double x :
       std::vector<double>{-112.8, -15, 0, 0.5, 1.0, 2.0, 3.9, 157.2}) {
    for (double y : std::vector<double>{0.5, 1.0, 2.0, 3.9, 912.37}) {
      stan::test::expect_ad(f, x, y);
    }
  }

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
