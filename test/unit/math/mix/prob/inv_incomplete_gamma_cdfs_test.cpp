#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

// Each function at two points: z > a + 1 (continued fraction) and
// z < a + 1 (power series side) of the incomplete gamma function.

TEST_F(AgradRev, mathMixScalFun_inv_gamma_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f, 0.4, 2.0, 1.5);
    stan::test::expect_ad(f, 2.0, 3.0, 1.0);
  };
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::inv_gamma_cdf(y, a, b);
  });
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::inv_gamma_lcdf(y, a, b);
  });
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::inv_gamma_lccdf(y, a, b);
  });
}

TEST_F(AgradRev, mathMixScalFun_inv_chi_square_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f, 0.1, 3.0);
    stan::test::expect_ad(f, 1.0, 4.0);
  };
  check([](const auto& y, const auto& nu) {
    return stan::math::inv_chi_square_cdf(y, nu);
  });
  check([](const auto& y, const auto& nu) {
    return stan::math::inv_chi_square_lcdf(y, nu);
  });
  check([](const auto& y, const auto& nu) {
    return stan::math::inv_chi_square_lccdf(y, nu);
  });
}

TEST_F(AgradRev, mathMixScalFun_scaled_inv_chi_square_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f, 0.5, 3.0, 1.0);
    stan::test::expect_ad(f, 2.0, 4.0, 0.8);
  };
  check([](const auto& y, const auto& nu, const auto& s) {
    return stan::math::scaled_inv_chi_square_cdf(y, nu, s);
  });
  check([](const auto& y, const auto& nu, const auto& s) {
    return stan::math::scaled_inv_chi_square_lcdf(y, nu, s);
  });
  check([](const auto& y, const auto& nu, const auto& s) {
    return stan::math::scaled_inv_chi_square_lccdf(y, nu, s);
  });
}
