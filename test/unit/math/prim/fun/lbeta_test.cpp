#include <test/unit/math/expect_near_rel.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, lbeta) {
  using stan::math::lbeta;

  EXPECT_FLOAT_EQ(0.0, lbeta(1.0, 1.0));
  EXPECT_FLOAT_EQ(2.981361, lbeta(0.1, 0.1));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(3.0, 4.0));
  EXPECT_FLOAT_EQ(-4.094345, lbeta(4.0, 3.0));
}

TEST(MathFunctions, lbeta_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, 1.0)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(1.0, nan)));

  EXPECT_TRUE(std::isnan(stan::math::lbeta(nan, nan)));
}

namespace lbeta_test_internal {
struct TestValue {
  double x;
  double y;
  double val;
};

// Test values generated in Mathematice, reproducible notebook at
// https://www.wolframcloud.com/obj/martin.modrak/Published/lbeta.nb
// Mathematica Code reproduced below for convenience:
//
// lbeta[x_,y_]:= LogGamma[x] + LogGamma[y] - LogGamma[x + y]
// out = OpenWrite["lbeta_test.txt"]
// xs= {8*10^-8,4*10^-3,1,1+10^-8,5,23, 19845};
// ys =  {7*10^-11,2*10^-5,1+10^-12,1/2,2,1624};
//  WriteString[out, "std::vector<TestValue> testValues = {"];
//  For[i = 1, i <= Length[xs], i++, {
//        For[j = 1, j <= Length[ys], j++, {
//      cx = xs[[i]];
//      cy = ys[[j]];

//     val = N[lbeta[cx,cy],24];
//          WriteString[out,"  {",CForm[cx],",",CForm[cy],",",
//            CForm[val],"},"]
// }]
// }]
// extremeXs = {3*10^15,10^20};
// lowYs = {3, 100, 12895};
//      For[i = 1, i <= Length[extremeXs], i++, {
//        For[j = 1, j <= Length[lowYs], j++, {
//      cx = extremeXs[[i]];
//          cy = lowYs[[j]];
//     val = N[lbeta[cx,cy],24];
//          WriteString[out,"  {",CForm[cx],".0,",CForm[cy],",",
//            CForm[val],"},"]
//      }]
//   }]
// WriteString[out,"};"];
//  Close[out];
//  FilePrint[%]
std::vector<TestValue> testValues = {
    {8.e-8, 7.e-11, 23.3834004912898500586445},
    {8.e-8, 0.00002, 16.3452312235394351410033},
    {8.e-8, 1.000000000001, 16.3412392022725295437606},
    {8.e-8, 0.5, 16.3412393131760679059067},
    {8.e-8, 2, 16.3412391222725327438921},
    {8.e-8, 1624, 16.3412385647081130254943},
    {0.004, 7.e-11, 23.3825258913787298259023},
    {0.004, 0.00002, 10.8247656947117878792194},
    {0.004, 1.000000000001, 5.52146091786223987264715},
    {0.004, 0.5, 5.52697992926150653113797},
    {0.004, 2, 5.51746889659270898022044},
    {0.004, 1624, 5.48959582574332555214719},
    {1, 7.e-11, 23.3825258738791892190926},
    {1, 0.00002, 10.8197782844102831106727},
    {1, 1.000000000001, -9.999999999995e-13},
    {1, 0.5, 0.693147180559945309417232},
    {1, 2, -0.693147180559945309417232},
    {1, 1624, -7.39264752072162326054032},
    {1.00000001, 7.e-11, 23.3825258738791892179411},
    {1.00000001, 0.00002, 10.8197782844099541286699},
    {1.00000001, 1.000000000001, -1.00009999500064491739816e-8},
    {1.00000001, 0.5, 0.693147174422888956122731},
    {1.00000001, 2, -0.693147195559945246917232},
    {1.00000001, 1624, -7.39264760042333353631934},
    {5, 7.e-11, 23.3825258737333558857627},
    {5, 0.00002, 10.8197366180283355258393},
    {5, 1.000000000001, -1.60943791243638370793409},
    {5, 0.5, -0.207395194346070587158746},
    {5, 2, -3.40119738166215537541324},
    {5, 1624, -33.7913357290267948074624},
    {23, 7.e-11, 23.3825258736208322915813},
    {23, 0.00002, 10.819704468465374949026},
    {23, 1.000000000001, -3.13549421593288398231784},
    {23, 0.5, -0.989947810259228199543883},
    {23, 2, -6.31354804627709531045369},
    {23, 1624, -121.714785277510463870251},
    {19845, 7.e-11, 23.3825258731460863706715},
    {19845, 0.00002, 10.8195688267825637640878},
    {19845, 1.000000000001, -9.89570736522763869861762},
    {19845, 0.5, -4.37548244086806082919414},
    {19845, 2, -19.7914651196913525680177},
    {19845, 1624, -5756.4146766727238501215},
    {3000000000000000.0, 3, -106.219018870176440545578},
    {3000000000000000.0, 100, -3204.60466298830574639047},
    {3000000000000000.0, 12895, -350396.988955562106921852},
    {100000000000000000000.0, 3, -137.461958399082795731692},
    {100000000000000000000.0, 100, -4246.03598061851596930944},
    {100000000000000000000.0, 12895, -484689.557363950217404711},
};
}  // namespace lbeta_test_internal

TEST(MathFunctions, lbeta_precomputed) {
  using lbeta_test_internal::TestValue;
  using lbeta_test_internal::testValues;
  using stan::test::expect_near_rel;

  for (TestValue t : testValues) {
    std::ostringstream msg;
    msg << std::setprecision(22) << "x = " << t.x << ", y = " << t.y;

    double val = stan::math::lbeta(t.x, t.y);
    expect_near_rel(msg.str(), val, t.val);
  }
}
