#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST_F(AgradRev, mathMixScalFun_von_mises_cdf) {
  auto f = [](const auto& x, const auto& mu, const auto& k) {
    return stan::math::von_mises_cdf(x, mu, k);
  };

  // large k
  stan::test::expect_ad(f, 0.9, 0.9, 100);

  // double instantiations
  stan::test::expect_ad(f, 10.6, 10.3, 0.5);
  stan::test::expect_ad(f, 0.9, 0.0, 0.9);
  stan::test::expect_ad(f, 3.0, 0.0, 0.4);
  stan::test::expect_ad(f, -3.0, 0.0, 0.4);

  // integer instantiations
  stan::test::expect_ad(f, 3, 0, 0.2);
  stan::test::expect_ad(f, -3, 2.0, 0.2);
  stan::test::expect_ad(f, -3.0, 0, 0.2);
  stan::test::expect_ad(f, 3.0, 2.0, 0.2);

  // nan instantiations
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 0.6, 0.3, nan);
  stan::test::expect_ad(f, 0.6, nan, 0.5);
  stan::test::expect_ad(f, 0.6, nan, nan);
  stan::test::expect_ad(f, nan, 0.3, 0.5);
  stan::test::expect_ad(f, nan, 0.3, nan);
  stan::test::expect_ad(f, nan, nan, 0.5);
  stan::test::expect_ad(f, nan, nan, nan);

  // check helper functions
  auto f2 = [](const auto& x, const auto& k) {
    return stan::math::internal::von_mises_cdf_centered(x, k);
  };
  stan::test::expect_ad(f2, 0.1, 12.7);

  auto f3 = [](const auto& x, const auto& k) {
    return stan::math::internal::von_mises_cdf_normalapprox(x, k);
  };
  stan::test::expect_ad(f3, 0.7, 62.7);

  auto f4 = [](const auto& x, const auto& k) {
    return stan::math::internal::von_mises_cdf_series(x, k);
  };
  stan::test::expect_ad(f4, -1.7, 12.3);
}
