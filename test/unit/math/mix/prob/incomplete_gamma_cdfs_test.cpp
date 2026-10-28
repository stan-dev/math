#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

// Each function at a point with z > a + 1 (continued fraction) and a
// point with z < a + 1 (power series side) of the incomplete gamma
// function; gamma and chi_square also at a point next to z = a + 1.

TEST_F(AgradRev, mathMixScalFun_gamma_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f, 4.0, 2.0, 1.5);
    stan::test::expect_ad(f, 1.0, 4.0, 1.5);
    stan::test::expect_ad(f, 2.5, 3.0, 1.5);
  };
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::gamma_cdf(y, a, b);
  });
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::gamma_lcdf(y, a, b);
  });
  check([](const auto& y, const auto& a, const auto& b) {
    return stan::math::gamma_lccdf(y, a, b);
  });
}

TEST_F(AgradRev, mathMixScalFun_chi_square_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f, 8.0, 3.0);
    stan::test::expect_ad(f, 2.0, 6.0);
    stan::test::expect_ad(f, 5.0, 3.0);
  };
  check([](const auto& y, const auto& nu) {
    return stan::math::chi_square_cdf(y, nu);
  });
  check([](const auto& y, const auto& nu) {
    return stan::math::chi_square_lcdf(y, nu);
  });
  check([](const auto& y, const auto& nu) {
    return stan::math::chi_square_lccdf(y, nu);
  });
}

TEST_F(AgradRev, mathMixScalFun_poisson_cdfs) {
  // the count is bound into the functor; only lambda is differentiated
  auto check = [](const auto& f) {
    stan::test::expect_ad(f(3), 6.5);
    stan::test::expect_ad(f(5), 2.0);
  };
  check([](int n) {
    return
        [n](const auto& lambda) { return stan::math::poisson_cdf(n, lambda); };
  });
  check([](int n) {
    return
        [n](const auto& lambda) { return stan::math::poisson_lcdf(n, lambda); };
  });
  check([](int n) {
    return [n](const auto& lambda) {
      return stan::math::poisson_lccdf(n, lambda);
    };
  });
}
