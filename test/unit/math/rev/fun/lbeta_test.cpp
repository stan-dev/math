#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev/core/nested_rev_autodiff.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <functional>

namespace lbeta_test_internal {
// TODO(martinmodrak) the function here should be replaced by helpers for
// testing identities once those are available

struct identity_tolerances {
  stan::test::relative_tolerance value;
  stan::test::relative_tolerance gradient;
};

template <class F1, class F2>
void expect_identity(const std::string& msg,
                     const identity_tolerances& tolerances, const F1 lh,
                     const F2 rh, double x_dbl, double y_dbl) {
  using stan::math::var;
  using stan::test::expect_near_rel;

  stan::math::nested_rev_autodiff nested;

  var x(x_dbl);
  var y(y_dbl);

  std::vector<var> vars = {x, y};

  var left = lh(x, y);
  double left_dbl = value_of(left);
  std::vector<double> gradients_left;
  left.grad(vars, gradients_left);

  nested.set_zero_all_adjoints();

  var right = rh(x, y);
  double right_dbl = value_of(right);
  std::vector<double> gradients_right;
  right.grad(vars, gradients_right);

  std::stringstream args;
  args << std::setprecision(22) << "args = [" << x << "," << y << "]";
  expect_near_rel(std::string() + args.str() + std::string(": ") + msg,
                  left_dbl, right_dbl, tolerances.value);

  for (size_t i = 0; i < gradients_left.size(); ++i) {
    std::stringstream grad_msg;
    grad_msg << "grad_" << i << ", " << args.str() << ": " << msg;
    expect_near_rel(grad_msg.str(), gradients_left[i], gradients_right[i],
                    tolerances.gradient);
  }
}
}  // namespace lbeta_test_internal

TEST(MathFunctions, lbeta_identities_gradient) {
  using stan::math::lbeta;
  using stan::math::pi;
  using stan::math::var;
  using stan::test::expect_near_rel;

  std::vector<double> to_test
      = {1e-100, 1e-8, 1e-1, 1, 2, 1 + 1e-6, 1e3, 1e30, 1e100};

  lbeta_test_internal::identity_tolerances tol{{1e-15, 1e-15}, {1e-10, 1e-10}};

  // All identities from https://en.wikipedia.org/wiki/Beta_function#Properties
  // Successors: beta(a,b) = beta(a + 1, b) + beta(a, b + 1)
  for (double x : to_test) {
    for (double y : to_test) {
      auto rh = [](const var& a, const var& b) {
        return stan::math::log_sum_exp(lbeta(a + 1, b), lbeta(a, b + 1));
      };
      lbeta_test_internal::expect_identity(
          "succesors", tol, static_cast<var (*)(const var&, const var&)>(lbeta),
          rh, x, y);
    }
  }

  // Sin: beta(x, 1 - x) == pi / sin(pi * x)
  for (double x : to_test) {
    if (x < 1) {
      std::stringstream msg;
      msg << std::setprecision(22) << "sin: x = " << x;
      double lh = lbeta(x, 1.0 - x);
      double rh = log(pi()) - log(sin(pi() * x));
      expect_near_rel(msg.str(), lh, rh, tol.value);
    }
  }

  // Inv: beta(1, x) == 1 / x
  for (double x : to_test) {
    std::stringstream msg;
    msg << std::setprecision(22) << "inv: x = " << x;
    double lh = lbeta(x, 1.0);
    double rh = -log(x);
    expect_near_rel(msg.str(), lh, rh, tol.value);
  }
}

namespace lbeta_test_internal {
struct TestValue {
  double x;
  double y;
  double val;
  double dx;
  double dy;
};

const double NaN = std::numeric_limits<double>::quiet_NaN();
// Test values generated in Mathematica, reproducible notebook at
// https://www.wolframcloud.com/obj/martin.modrak/Published/lbeta.nb
// Mathematica Code reproduced below for convenience:
//
// toCString[x_] := ToString[CForm[N[x, 24]]];
// singleTest[x_,y_]:= Module[{val, cdx,cdy, big},{
// val = toCString[lbeta[x,y]];
// cdx = If[x > 10^6 || y > 10^6,"NaN", toCString[dlbetadx[x,y]]];
// cdy = If[x > 10^6 || y > 10^6,"NaN", toCString[dlbetady[x,y]]];
// StringJoin["  {",toCString[x],",",toCString[y],",",
//            val,",",cdx,",",cdy,"},\n"]
// }];
// xs= SetPrecision[{0.00000001,0.004,1,1.00000001,5.0,23.0, 19845.0},
//   Infinity];
// ys = SetPrecision[{0.00000000001,0.00002,1.000000000001,0.5,2.0,1624.0},
//   Infinity];
// out = "std::vector<TestValue> testValues = {\n";
//  For[i = 1, i <= Length[xs], i++, {
//        For[j = 1, j <= Length[ys], j++, {
//      out = StringJoin[out, singleTest[xs[[i]],ys[[j]] ]];
// }]
// }]
// extremeXs = {3*10^15+10^-1,10^20 + 10^-1};
// lowYs = {3, 100, 12895};
//      For[i = 1, i <= Length[extremeXs], i++, {
//        For[j = 1, j <= Length[lowYs], j++, {
//     out = StringJoin[out, singleTest[extremeXs[[i]],lowYs[[j]] ]];
//      }]
//   }]
// out = StringJoin[out,"};\n"];
// out

std::vector<TestValue> testValues = {
    {1.00000000000000002092256e-8, 9.9999999999999993949697e-12,
     25.3294355232675861176219, -99900.0999000999063327762,
     -9.99000999000999061687346e10},
    {1.00000000000000002092256e-8, 0.0000200000000000000016360611,
     18.4211806189936875171022, -9.99500249875391429884343e7,
     -24.987506263325417893104},
    {1.00000000000000002092256e-8, 1.00000000000100008890058,
     18.4206807439523654512049, -9.9999999999999997909389e7,
     -1.64493405482525322002407e-8},
    {1.00000000000000002092256e-8, 0.5, 18.4206807578153088979269,
     -9.99999986137056696865339e7, -4.93480211640069781497916e-8},
    {1.00000000000000002092256e-8, 2., 18.4206807339523655012214,
     -1.00000000999999987907744e8, -6.44934064827657426602468e-9},
    {1.00000000000000002092256e-8, 1624., 18.4206806642568128109103,
     -1.00000007969555253717259e8, -6.15953168081550363455631e-12},
    {0.00400000000000000008326673, 9.9999999999999993949697e-12,
     25.3284360254344369758523, -6.25016352130489628103017e-7,
     -9.99999997500065672476642e10},
    {0.00400000000000000008326673, 0.0000200000000000000016360611,
     10.8247656947117877977197, -1.24381380143768420568693,
     -49751.2503414755942898455},
    {0.00400000000000000008326673, 1.00000000000100008890058,
     5.52146091786223985124723, -250.000000000001630310514,
     -0.00656057236120388443368478},
    {0.00400000000000000008326673, 0.5, 5.52697992926150651043563,
     -248.626750675814271729954, -0.0196056092953810030104119},
    {0.00400000000000000008326673, 2., 5.51746889659270895932082,
     -250.996015936254974875428, -0.00257650861619352965053606},
    {0.00400000000000000008326673, 1624., 5.48959582574332553066745,
     -257.962997163701114672183, -2.46380963715242861602996e-6},
    {1., 9.9999999999999993949697e-12, 25.3284360229345025847009,
     -1.64493406683620576791743e-11, -1.00000000000000006050303e11},
    {1., 0.0000200000000000000016360611, 10.8197782844102830288697,
     -0.0000328982005228616875473718, -49999.9999999999959098473},
    {1., 1.00000000000100008890058, -1.00008890058184092270907e-12,
     -1.00000000000064499140186, -0.99999999999899991109942},
    {1., 0.5, 0.693147180559945309417232, -0.613705638880109381165536, -2.},
    {1., 2., -0.693147180559945309417232, -1.5, -0.5},
    {1., 1624., -7.39264752072162326054032, -7.97017103579949420353338,
     -0.000615763546798029556650246},
    {1.00000000999999993922529, 9.9999999999999993949697e-12,
     25.3284360229345025845364, -1.64493404279506817585652e-11,
     -1.00000000000000006066752e11},
    {1.00000000999999993922529, 0.0000200000000000000016360611,
     10.8197782844099540468688, -0.0000328982000420519232463277,
     -50000.0000000164447694858},
    {1.00000000999999993922529, 1.00000000000100008890058,
     -1.00009999781323229867408e-8, -0.999999990000645152172528,
     -1.00000000644834052017627},
    {1.00000000999999993922529, 0.5, 0.693147174422888993420513,
     -0.613705631778790840053979, -2.00000000934802190719463},
    {1.00000000999999993922529, 2., -0.693147195559945155755169,
     -1.49999998750000018846838, -0.500000003949340636774571},
    {1.00000000999999993922529, 1624., -7.39264760042333305193451,
     -7.97017101935630949522315, -0.000615763552953769552598789},
    {5., 9.9999999999999993949697e-12, 25.3284360229136692513677,
     -2.21322955736871363309369e-12, -1.00000000002083339383622e11},
    {5., 0.0000200000000000000016360611, 10.8197366180283354440329,
     -4.42644935682442819025907e-6, -50002.0833048615780771434},
    {5., 1.00000000000100008890058, -1.60943791243638391092375,
     -0.200000000000181339075453, -2.28333333333186959210634},
    {5., 0.5, -0.207395194346070587158746, -0.104975480149950651006806,
     -3.5746031746031746031746},
    {5., 2., -3.40119738166215537541324, -0.366666666666666666666667, -1.45},
    {5., 1624., -33.7913357290267948074624, -5.88929697199582686904279,
     -0.00307503307646402839938757},
    {23., 9.9999999999999993949697e-12, 25.3284360228975944521988,
     -4.44371335365687819119143e-13, -1.00000000003690819300504e11},
    {23., 0.0000200000000000000016360611, 10.8197044684653748672169,
     -8.88742275864946586539857e-7, -50003.6907812407549378356},
    {23., 1.00000000000100008890058, -3.13549421593288431429853,
     -0.0434782608696077679481057, -3.7342915110852376727225},
    {23., 0.5, -0.989947810259228199543883, -0.0219753695482036950492042,
     -5.09908298088536929896051},
    {23., 2., -6.31354804627709531045369, -0.085144927536231884057971,
     -2.77595817775350686913481},
    {23., 1624., -121.714785277510463870251, -4.29280953187037223167221,
     -0.0140675098349510427723158},
    {19845., 9.9999999999999993949697e-12, 25.3284360228297736063546,
     -5.03917962049126253574903e-16, -1.0000000001047290388493e11},
    {19845., 0.0000200000000000000016360611, 10.8195688267825636822676,
     -1.00783592359038628000151e-9, -50010.4728649374470305881},
    {19845., 1.00000000000100008890058, -9.89570736522763962966881,
     -0.000050390526581053165215596, -10.4729482251687436517705},
    {19845., 0.5, -4.37548244086806082919414, -0.0000251955806911474238396687,
     -11.859217391344389445166},
    {19845., 2., -19.7914651196913525680177, -0.000100778514084781870541745,
     -9.47299861315789246077826},
    {19845., 1624., -5756.4146766727238501215, -0.078659853852481428671324,
     -2.58200241624359293360231},
    {3.0000000000000001e15, 3., -106.219018870176440645578, NaN, NaN},
    {3.0000000000000001e15, 100., -3204.6046629883057497238, NaN, NaN},
    {3.0000000000000001e15, 12895., -350396.988955562107351686, NaN, NaN},
    {1.000000000000000000001e20, 3., -137.461958399082795731695, NaN, NaN},
    {1.000000000000000000001e20, 100., -4246.03598061851596930954, NaN, NaN},
    {1.000000000000000000001e20, 12895., -484689.557363950217404724, NaN, NaN},
};
}  // namespace lbeta_test_internal

TEST(MathFunctions, lbeta_precomputed) {
  using lbeta_test_internal::TestValue;
  using lbeta_test_internal::testValues;
  using stan::math::is_nan;
  using stan::math::value_of;
  using stan::math::var;

  for (TestValue t : testValues) {
    std::ostringstream msg;
    msg << std::setprecision(22) << "x = " << t.x << ", y = " << t.y;

    var x(t.x);
    var y(t.y);
    var val = stan::math::lbeta(x, y);

    std::vector<var> vars;
    vars.push_back(x);
    vars.push_back(y);

    std::vector<double> gradients;
    val.grad(vars, gradients);

    for (int i = 0; i < 2; ++i) {
      EXPECT_FALSE(is_nan(gradients[i]));
    }

    double tol = std::max(fabs(t.val) * 1e-15, 1e-16);
    EXPECT_NEAR(value_of(val), t.val, tol) << msg.str();

    std::function<double(double)> tol_grad;
    if (x < 1e-4 || y < 1e-4) {
      tol_grad = [](double x) { return std::max(fabs(x) * 1e-8, 1e-7); };
    } else {
      tol_grad = [](double x) { return std::max(fabs(x) * 1e-10, 1e-8); };
    }
    if (!is_nan(t.dx)) {
      EXPECT_NEAR(gradients[0], t.dx, tol_grad(t.dx)) << "dx: " << msg.str();
    }
    if (!is_nan(t.dy)) {
      EXPECT_NEAR(gradients[1], t.dy, tol_grad(t.dy)) << "dy: " << msg.str();
    }
  }
}
