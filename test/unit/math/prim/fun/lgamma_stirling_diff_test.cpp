#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <cmath>

TEST(MathFunctions, lgamma_stirling_diff_errors_special_cases) {
  using stan::math::lgamma_stirling_diff;

  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();

  EXPECT_TRUE(std::isnan(lgamma_stirling_diff(nan)));
  EXPECT_FLOAT_EQ(lgamma_stirling_diff(inf), 0);
  EXPECT_THROW(std::isnan(lgamma_stirling_diff(-1.0)), std::domain_error);
  EXPECT_FLOAT_EQ(lgamma_stirling_diff(0.0), inf);
}

TEST(MathFunctions, lgamma_stirling_diff_accuracy) {
  using stan::math::lgamma_stirling_diff;
  using stan::math::lgamma_stirling_diff_useful;

  double start = std::nextafter(10, 11);
  for (double x = start; x < 1e8; x *= 1.5) {
    long double x_l = static_cast<long double>(x);
    long double stirling
        = 0.5 * std::log(2 * static_cast<long double>(stan::math::pi()))
          + (x_l - 0.5) * std::log(x_l) - x_l;
    long double lgamma_res = std::lgamma(x_l);
    double diff_actual = static_cast<double>(lgamma_res - stirling);
    double diff = lgamma_stirling_diff(x);
    double tol = std::max(1e-06 * 0.5 * (fabs(diff_actual) + fabs(diff)), 1e-6);

    EXPECT_NEAR(diff, diff_actual, tol)  << "x = " << x 
      << "; lgamma = " << lgamma_res << "; stirling = " << stirling;
  }

  double before_useful = std::nextafter(lgamma_stirling_diff_useful, 0);
  double after_useful = std::nextafter(lgamma_stirling_diff_useful,
                                    stan::math::positive_infinity());

  double value_before = lgamma_stirling_diff(before_useful);                                  
  double value_at = lgamma_stirling_diff(lgamma_stirling_diff_useful);
  double value_after = lgamma_stirling_diff(after_useful);
  double diff_before = value_at - value_before;
  double diff_after = value_after - value_at;
  double tol = std::max(
      1e-14 * (0.5 * (fabs(diff_before) + fabs(diff_after))), 1e-14);
  EXPECT_NEAR(diff_before, diff_after, tol) << "diff around useful cutoff";
}

namespace lgamma_stirling_diff_test_internal {
struct TestValue {
  double x;
  double val;
};

std::vector<TestValue> testValues = {
    {1.049787068367863943, 0.077388806767834476832},
    {1.1353352832366126919, 0.071790358566585005482},
    {1.3678794411714423216, 0.059960812482712981438},
    {2., 0.041340695955409294094},
    {3.7182818284590452354, 0.022358812123082674471},
    {8.3890560989306502272, 0.0099288907523535997267},
    {21.085536923187667741, 0.0039518599801395734578},
    {55.598150033144239078, 0.0014988346688724404687},
    {149.41315910257660342, 0.0005577367442531155476},
    {404.42879349273512261, 0.00020605188772717995062},
    {1097.6331584284585993, 0.000075920930766205666598},
    {2981.9579870417282747, 0.00002794584410078046085},
    {8104.0839275753840077, 0.000010282881326966996581},
    {22027.465794806716517, 3.7831557249429676373e-6},
    {59875.141715197818455, 1.3917851539949910276e-6},
    {162755.79141900392081, 5.1201455018391551878e-7},
    {442414.39200892050333, 1.8836035815859912686e-7},
    {1.2026052841647767777e6, 6.9294002305342748064e-8},
    {3.2690183724721106393e6, 2.5491852243801980915e-8},
    {8.8861115205078726368e6, 9.3779301712579119232e-9},
    {2.4154953753575298215e7, 3.4499479561619419703e-9},
    {6.5659970137330511139e7, 1.2691649593966957356e-9},
    {1.7848230196318726084e8, 4.6689970051216164913e-10},
    {4.8516519640979027797e8, 1.7176280151585029843e-10},
    {1.3188157354832146972e9, 6.3188003518019870566e-11},
    {3.5849128471315915617e9, 2.3245567434090096532e-11},
    {9.7448034472489026e9, 8.5515663588740238069e-12},
    {2.6489122130843472294e10, 3.1459454534471511195e-12},
    {7.2004899338385872524e10, 1.1573286553975954675e-12},
    {1.9572960942983876427e11, 4.2575741900310182743e-13},
    {5.3204824060279861668e11, 1.5662740137796255552e-13},
    {1.4462570642924751737e12, 5.7620000891128517638e-14},
};

}  // namespace lgamma_stirling_diff_test_internal

TEST(MathFunctions, lgamma_stirling_diff_precomputed) {
  using lgamma_stirling_diff_test_internal::TestValue;
  using lgamma_stirling_diff_test_internal::testValues;
  using stan::math::lgamma_stirling_diff;
  using stan::math::lgamma_stirling_diff_useful;

  for (TestValue t : testValues) {
    double tol = std::max(1e-16 * fabs(t.val), 1e-17);
    if (t.x < lgamma_stirling_diff_useful) {
      tol = std::max(1e-14 * fabs(t.val), 1e-14);
    }
    EXPECT_NEAR(lgamma_stirling_diff(t.x), t.val, tol) << "x = " << t.x;
  }
}
