#include <test/unit/math/test_ad.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

TEST(mathMixScalFun, lbeta) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lbeta;
    return lbeta(x1, x2);
  };

  stan::test::expect_common_nonzero_binary(f);

  stan::test::expect_ad(f, 0.5, 1.2);
  stan::test::expect_ad(f, 3.0, 6.0);
  stan::test::expect_ad(f, 3.1, 2.5);
  stan::test::expect_ad(f, 5.2, 6.7);
  stan::test::expect_ad(f, 6.7, 3.1);
  stan::test::expect_ad(f, 7.5, 1.8);
  stan::test::expect_ad(f, 12.3, 4.8);
}

TEST(mathMixScalFun, lbeta_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::lbeta;
    return lbeta(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

// The partial of lbeta(a, b) in an argument x is digamma(x) - digamma(a + b).
// For a large argument the two digamma values agree in most digits, and the
// plain difference keeps no correct digits of x * d/dx. References from
// mpmath at 120 digits (closed-form digamma differences, checked against
// mp.diff). The derivatives are compared as x * d/dx, in both argument
// orders, for var and fvar<double>.
TEST_F(AgradRev, mathMixScalFun_lbeta_large_argument_derivatives) {
  using stan::math::fvar;
  using stan::math::var;
  struct TestValue {
    double a;
    double b;
    double value;
    double d_a;
    double d_b;
  };
  const std::vector<TestValue> test_values = {
      {0x1.0000000000000p-1, 0x1.6bcc41e900000p+46, -15.545730708033618,
       -34.199701327938065, -5.0000000000000126e-15},
      {0x1.8000000000000p+0, 0x1.6bcc41e900000p+46, -48.475069190510212,
       -32.199701327938072, -1.4999999999999961e-14},
      {0x1.0000000000000p-1, 0x1.7d78400000000p+26, -8.6379754278014822,
       -20.384190769973788, -5.0000000125000003e-09},
      {0x1.9000000000000p+3, 0x1.2a05f20000000p+33, -269.08878911950677,
       -20.540655279865543, -1.24999999928125e-09},
  };
  auto expect_log_scale
      = [](double x, double got, double expected, const std::string& where) {
          const double scaled = x * expected;
          EXPECT_NEAR(x * got, scaled, 1e-11 * std::max(1.0, std::fabs(scaled)))
              << where;
        };
  for (const auto& t : test_values) {
    for (bool swap : {false, true}) {
      const double x1 = swap ? t.b : t.a;
      const double x2 = swap ? t.a : t.b;
      const double d1 = swap ? t.d_b : t.d_a;
      const double d2 = swap ? t.d_a : t.d_b;
      const std::string where
          = "lbeta(" + std::to_string(x1) + ", " + std::to_string(x2) + ")";

      var v1 = x1;
      var v2 = x2;
      var lp = stan::math::lbeta(v1, v2);
      lp.grad();
      EXPECT_NEAR(lp.val(), t.value, 1e-13 * std::fabs(t.value)) << where;
      expect_log_scale(x1, v1.adj(), d1, "var d/dx1, " + where);
      expect_log_scale(x2, v2.adj(), d2, "var d/dx2, " + where);
      stan::math::recover_memory();

      fvar<double> f1(x1, 1.0);
      fvar<double> f2(x2, 0.0);
      expect_log_scale(x1, stan::math::lbeta(f1, f2).d_, d1,
                       "fvar d/dx1, " + where);
      fvar<double> g1(x1, 0.0);
      fvar<double> g2(x2, 1.0);
      expect_log_scale(x2, stan::math::lbeta(g1, g2).d_, d2,
                       "fvar d/dx2, " + where);
    }
  }
}
