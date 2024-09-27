#include <test/unit/math/test_ad.hpp>
#include <stan/math/prim/core/complex_base.hpp>

#include <vector>

TEST(mixScalFun, norm) {
  auto f = [](const auto& x) { return norm(x); };
  stan::test::expect_complex_common(f);
}
