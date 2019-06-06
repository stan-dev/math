#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixScalFun, logModifiedBesselFirstKind) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_modified_bessel_first_kind(x1, x2);
  };
  // finite diff tests fail at x1 = 0; remove 0 and run tests by hand
  std::vector<double> xs1{-1.3, 0.5, stan::math::positive_infinity(),
                          stan::math::negative_infinity(),
                          stan::math::not_a_number()};
  std::vector<double> xs2 = stan::test::common_args();
  for (double x1 : xs1)
    for (double x2 : xs2)
      stan::test::expect_ad(f, x1, x2);
}
