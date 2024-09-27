#include <test/unit/math/test_ad.hpp>
#include <stan/math/prim/core/complex_base.hpp>

#include <vector>

TEST(mixScalFun, proj) {
  auto f = [](const auto& x) { return proj(x); };
  stan::test::expect_complex_common(f);
}
