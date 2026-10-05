#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(MathMixProbWienerLpdf, fiveParamWGradientExpectAd) {
  auto f = [](const auto& w) {
    return stan::math::wiener_lpdf(6.0, 10.0, 0.01, w, -3.0, 0.2);
  };

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-5;
  stan::test::expect_ad(tols, f, 0.1);
}

TEST(MathMixProbWienerLpdf, fiveParamZeroSvWGradientExpectAd) {
  auto f = [](const auto& w) {
    return stan::math::wiener_lpdf(6.0, 10.0, 0.01, w, -3.0, 0.0);
  };

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-5;
  stan::test::expect_ad(tols, f, 0.1);
}

TEST(MathMixProbWienerLpdf, fullParamWGradientExpectAd) {
  auto f = [](const auto& w) {
    return stan::math::wiener_lpdf(6.0, 10.0, 0.01, w, -3.0, 0.2, 0.1, 0.0);
  };

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-5;
  stan::test::expect_ad(tols, f, 0.1);
}

TEST(MathMixProbWienerLpdf, existingFullRowsWGradientExpectAd) {
  struct Case {
    const char* name;
    double y;
    double a;
    double t0;
    double w;
    double v;
    double sv;
    double sw;
    double st0;
  };

  const std::vector<Case> cases = {
      {"row_0", 2.0, 2.0, 1e-9, 0.10, 2.0, 0.0, 0.00, 0.000},
      {"row_1", 3.0, 2.0, 0.01, 0.50, 2.0, 0.2, 0.00, 0.000},
      {"row_2", 4.0, 10.0, 0.01, 0.80, 4.0, 0.0, 0.10, 0.000},
      {"row_3", 5.0, 4.0, 0.01, 0.70, 3.0, 0.0, 0.00, 0.007},
      {"row_4", 6.0, 10.0, 0.01, 0.10, -3.0, 0.2, 0.10, 0.000},
      {"row_5", 7.0, 1.0, 0.01, 0.90, 1.0, 0.2, 0.00, 0.007},
      {"row_6", 8.0, 3.0, 0.01, 0.70, -1.0, 0.0, 0.10, 0.007},
      {"row_7", 8.85, 1.7, 0.01, 0.92, -7.3, 0.7, 0.01, 0.009},
      {"row_8", 8.9, 2.4, 0.01, 0.90, -4.9, 0.0, 0.00, 0.009},
      {"row_9", 9.0, 11.0, 0.01, 0.12, 4.5, 0.7, 0.10, 0.009},
      {"row_10", 1.0, 1.5, 0.10, 0.50, 3.0, 0.5, 0.20, 0.000},
  };

  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-4;

  for (const auto& c : cases) {
    SCOPED_TRACE(c.name);

    auto f = [c](const auto& w) {
      return stan::math::wiener_lpdf(c.y, c.a, c.t0, w, c.v, c.sv, c.sw, c.st0);
    };

    // The row sweep is intended to check the reverse-mode w adjoint against
    // finite differences. Some full-Wiener rows are not stable enough for
    // higher-order mixed-mode finite-difference checks
    stan::test::expect_ad<true>(tols, f, c.w);
  }
}
