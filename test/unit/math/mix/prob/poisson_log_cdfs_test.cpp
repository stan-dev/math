#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

// Each function at two points: exp(alpha) > n + 2 (continued fraction) and
// exp(alpha) < n + 2 (power series side) of the incomplete gamma function.
// The count is bound into the functor; only alpha is differentiated.
TEST_F(AgradRev, mathMixScalFun_poisson_log_cdfs) {
  auto check = [](const auto& f) {
    stan::test::expect_ad(f(3), 2.0);
    stan::test::expect_ad(f(5), 0.5);
  };
  check([](int n) {
    return [n](const auto& alpha) {
      return stan::math::poisson_log_cdf(n, alpha);
    };
  });
  check([](int n) {
    return [n](const auto& alpha) {
      return stan::math::poisson_log_lcdf(n, alpha);
    };
  });
  check([](int n) {
    return [n](const auto& alpha) {
      return stan::math::poisson_log_lccdf(n, alpha);
    };
  });
}
