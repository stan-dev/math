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

    EXPECT_NEAR(diff, diff_actual, tol)
        << "x = " << x << "; lgamma = " << lgamma_res
        << "; stirling = " << stirling;
  }

  double before_useful = std::nextafter(lgamma_stirling_diff_useful, 0);
  double after_useful = std::nextafter(lgamma_stirling_diff_useful,
                                       stan::math::positive_infinity());

  double value_before = lgamma_stirling_diff(before_useful);
  double value_at = lgamma_stirling_diff(lgamma_stirling_diff_useful);
  double value_after = lgamma_stirling_diff(after_useful);
  double diff_before = value_at - value_before;
  double diff_after = value_after - value_at;
  double tol
      = std::max(1e-14 * (0.5 * (fabs(diff_before) + fabs(diff_after))), 1e-14);
  EXPECT_NEAR(diff_before, diff_after, tol) << "diff around useful cutoff";
}

namespace lgamma_stirling_diff_test_internal {
struct TestValue {
  double x;
  double val;
};

// Test values precomputed in Mathematica, reproducible notebook available at
// https://www.wolframcloud.com/obj/martin.modrak/Published/lgamma_stirling_diff.nb
// Code reproduced also below for convenience
//
// lgs[x_]:=1/2 * Log[2*Pi] + (x - 1/2)*Log[x] -x
// sdiff[x_]:= LogGamma[x] - lgs[x]
// toCString[x_] := ToString[CForm[N[x, 24]]];
// stride = 1.0;
// xsraw= Table[1.0 +Exp[i/stride],{i,-3.0 * stride,28.0 * stride}];
// xs = SetPrecision[xsraw, Infinity];
//  out = "std::vector<TestValue> testValues = {\n";
//  For[i = 1, i <= Length[xs], i++, {
//          out = StringJoin[out,"  {",toCString[xs[[i]]],",",
//          toCString[sdiff[xs[[i]]]],"},\n"]
// }]
// out = StringJoin[out, "};\n"];
// out
std::vector<TestValue> testValues = {
    {1.04978706836786384037907, 0.0773888067678344840629031},
    {1.13533528323661281334012, 0.071790358566584998101908},
    {1.36787944117144233402428, 0.059960812482712980909299},
    {2., 0.0413406959554092940938221},
    {3.7182818284590450907956, 0.0223588121230826753367302},
    {8.38905609893065040694182, 0.00992889075235359951416742},
    {21.0855369231876679236848, 0.00395185998013957342350829},
    {55.5981500331442362039525, 0.00149883466887244054615435},
    {149.413159102576599934764, 0.000557736744253115560609879},
    {404.428793492735110248759, 0.000206051887727179956913261},
    {1097.63315842845850056619, 0.0000759209307662056734249356},
    {2981.95798704172830184689, 0.0000279458441007804605962937},
    {8104.08392757538422301877, 0.0000102828813269669963078472},
    {22027.4657948067178949714, 3.78315572494296740059554e-6},
    {59875.14171519781666575, 1.39178515399499106920527e-6},
    {162755.791419003915507346, 5.12014550183915535452136e-7},
    {442414.392008920491207391, 1.88360358158599132020887e-7},
    {1.20260528416477679274976e6, 6.92940023053427471994727e-8},
    {3.26901837247211067005992e6, 2.54918522438019806747147e-8},
    {8.88611152050787210464478e6, 9.37793017125791248474171e-9},
    {2.41549537535752989351749e7, 3.4499479561619418674483e-9},
    {6.56599701373305097222328e7, 1.26916495939669576293556e-9},
    {1.78482301963187247514725e8, 4.66899700512161684000123e-10},
    {4.85165196409790277481079e8, 1.71762801515850298601024e-10},
    {1.31881573548321461677551e9, 6.31880035180198744193528e-11},
    {3.584912847131591796875e9, 2.32455674340900950067326e-11},
    {9.74480344724890327453613e9, 8.55156635887402321498741e-12},
    {2.64891221308434715270996e10, 3.1459454534471512105765e-12},
    {7.20048993383858795166016e10, 1.1573286553975953551001e-12},
    {1.95729609429838775634766e11, 4.25757419003101802709993e-13},
    {5.32048240602798645019531e11, 1.56627401377962547173856e-13},
    {1.44625706429247509765625e12, 5.76200008911285206664666e-14},
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
