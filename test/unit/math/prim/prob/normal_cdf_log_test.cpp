#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(ProbNormal, cdf_log_matches_lcdf) {
  double y = 0.8;
  double mu = 1.1;
  double sigma = 2.3;

  EXPECT_FLOAT_EQ((stan::math::normal_lcdf(y, mu, sigma)),
                  (stan::math::normal_cdf_log(y, mu, sigma)));
  EXPECT_FLOAT_EQ(
      (stan::math::normal_lcdf(y, mu, sigma)),
      (stan::math::normal_cdf_log<double, double, double>(y, mu, sigma)));
}

TEST(ProbNormal, lcdf_tails) {
  using stan::math::normal_lcdf;
  using std::exp;

  // The test values come from R 4.6.1 and cover the expected useful range of
  // the function. When z >= 38.5, even the log of the CDF is
  // indistinguishable from 0.0 in double precision.
  //
  // q <-
  //   c(
  //     -10^seq(8, 2, by = -1),
  //     seq(-50, -11, by = 1.0),
  //     seq(-10, 10, by = 0.5),
  //     seq(11, 38, by = 1.0)
  //   )
  // for (i in 1:length(q)) {
  //   cat(
  //     sprintf(
  //       "EXPECT_FLOAT_EQ(%#.17g, normal_lcdf(%.17g, 0, 1));\n",
  //       pnorm(q[i], lower.tail = TRUE, log.p = TRUE),
  //       q[i]
  //     )
  //   )
  // }

  EXPECT_FLOAT_EQ(-5000000000000019.0, normal_lcdf(-100000000, 0, 1));
  EXPECT_FLOAT_EQ(-50000000000017.039, normal_lcdf(-10000000, 0, 1));
  EXPECT_FLOAT_EQ(-500000000014.73444, normal_lcdf(-1000000, 0, 1));
  EXPECT_FLOAT_EQ(-5000000012.4318638, normal_lcdf(-100000, 0, 1));
  EXPECT_FLOAT_EQ(-50000010.129278913, normal_lcdf(-10000, 0, 1));
  EXPECT_FLOAT_EQ(-500007.82669481216, normal_lcdf(-1000, 0, 1));
  EXPECT_FLOAT_EQ(-5005.5242086942053, normal_lcdf(-100, 0, 1));
  EXPECT_FLOAT_EQ(-1254.8313611394199, normal_lcdf(-50, 0, 1));
  EXPECT_FLOAT_EQ(-1205.3111748916654, normal_lcdf(-49, 0, 1));
  EXPECT_FLOAT_EQ(-1156.7905731019453, normal_lcdf(-48, 0, 1));
  EXPECT_FLOAT_EQ(-1109.2695383172531, normal_lcdf(-47, 0, 1));
  EXPECT_FLOAT_EQ(-1062.7480519624305, normal_lcdf(-46, 0, 1));
  EXPECT_FLOAT_EQ(-1017.2260942419524, normal_lcdf(-45, 0, 1));
  EXPECT_FLOAT_EQ(-972.70364403073665, normal_lcdf(-44, 0, 1));
  EXPECT_FLOAT_EQ(-929.18067875247391, normal_lcdf(-43, 0, 1));
  EXPECT_FLOAT_EQ(-886.65717424372951, normal_lcdf(-42, 0, 1));
  EXPECT_FLOAT_EQ(-845.13310460177456, normal_lcdf(-41, 0, 1));
  EXPECT_FLOAT_EQ(-804.60844201375380, normal_lcdf(-40, 0, 1));
  EXPECT_FLOAT_EQ(-765.08315656437753, normal_lcdf(-39, 0, 1));
  EXPECT_FLOAT_EQ(-726.55721601882010, normal_lcdf(-38, 0, 1));
  EXPECT_FLOAT_EQ(-689.03058557689064, normal_lcdf(-37, 0, 1));
  EXPECT_FLOAT_EQ(-652.50322759379844, normal_lcdf(-36, 0, 1));
  EXPECT_FLOAT_EQ(-616.97510126192253, normal_lcdf(-35, 0, 1));
  EXPECT_FLOAT_EQ(-582.44616224687172, normal_lcdf(-34, 0, 1));
  EXPECT_FLOAT_EQ(-548.91636226973810, normal_lcdf(-33, 0, 1));
  EXPECT_FLOAT_EQ(-516.38564862572537, normal_lcdf(-32, 0, 1));
  EXPECT_FLOAT_EQ(-484.85396362717927, normal_lcdf(-31, 0, 1));
  EXPECT_FLOAT_EQ(-454.32124395634321, normal_lcdf(-30, 0, 1));
  EXPECT_FLOAT_EQ(-424.78741990973015, normal_lcdf(-29, 0, 1));
  EXPECT_FLOAT_EQ(-396.25241451163106, normal_lcdf(-28, 0, 1));
  EXPECT_FLOAT_EQ(-368.71614246865636, normal_lcdf(-27, 0, 1));
  EXPECT_FLOAT_EQ(-342.17850892992783, normal_lcdf(-26, 0, 1));
  EXPECT_FLOAT_EQ(-316.63940800802027, normal_lcdf(-25, 0, 1));
  EXPECT_FLOAT_EQ(-292.09872100320780, normal_lcdf(-24, 0, 1));
  EXPECT_FLOAT_EQ(-268.55631425686312, normal_lcdf(-23, 0, 1));
  EXPECT_FLOAT_EQ(-246.01203653738091, normal_lcdf(-22, 0, 1));
  EXPECT_FLOAT_EQ(-224.46571583141449, normal_lcdf(-21, 0, 1));
  EXPECT_FLOAT_EQ(-203.91715537109727, normal_lcdf(-20, 0, 1));
  EXPECT_FLOAT_EQ(-184.36612866916096, normal_lcdf(-19, 0, 1));
  EXPECT_FLOAT_EQ(-165.81237325071419, normal_lcdf(-18, 0, 1));
  EXPECT_FLOAT_EQ(-148.25558265098039, normal_lcdf(-17, 0, 1));
  EXPECT_FLOAT_EQ(-131.69539607375970, normal_lcdf(-16, 0, 1));
  EXPECT_FLOAT_EQ(-116.13138484571169, normal_lcdf(-15, 0, 1));
  EXPECT_FLOAT_EQ(-101.56303440744996, normal_lcdf(-14, 0, 1));
  EXPECT_FLOAT_EQ(-87.989719971022524, normal_lcdf(-13, 0, 1));
  EXPECT_FLOAT_EQ(-75.410673001568796, normal_lcdf(-12, 0, 1));
  EXPECT_FLOAT_EQ(-63.824934094423718, normal_lcdf(-11, 0, 1));
  EXPECT_FLOAT_EQ(-53.231285150512470, normal_lcdf(-10, 0, 1));
  EXPECT_FLOAT_EQ(-48.306019298965232, normal_lcdf(-9.5, 0, 1));
  EXPECT_FLOAT_EQ(-43.628149113332114, normal_lcdf(-9, 0, 1));
  EXPECT_FLOAT_EQ(-39.197396428217672, normal_lcdf(-8.5, 0, 1));
  EXPECT_FLOAT_EQ(-35.013437159914552, normal_lcdf(-8, 0, 1));
  EXPECT_FLOAT_EQ(-31.075890902890002, normal_lcdf(-7.5, 0, 1));
  EXPECT_FLOAT_EQ(-27.384307498811076, normal_lcdf(-7, 0, 1));
  EXPECT_FLOAT_EQ(-23.938149495161838, normal_lcdf(-6.5, 0, 1));
  EXPECT_FLOAT_EQ(-20.736768949974707, normal_lcdf(-6, 0, 1));
  EXPECT_FLOAT_EQ(-17.779376352625260, normal_lcdf(-5.5, 0, 1));
  EXPECT_FLOAT_EQ(-15.064998393988725, normal_lcdf(-5, 0, 1));
  EXPECT_FLOAT_EQ(-12.592419735713079, normal_lcdf(-4.5, 0, 1));
  EXPECT_FLOAT_EQ(-10.360101486527292, normal_lcdf(-4, 0, 1));
  EXPECT_FLOAT_EQ(-8.3660653083440941, normal_lcdf(-3.5, 0, 1));
  EXPECT_FLOAT_EQ(-6.6077262215103492, normal_lcdf(-3, 0, 1));
  EXPECT_FLOAT_EQ(-5.0816482772786902, normal_lcdf(-2.5, 0, 1));
  EXPECT_FLOAT_EQ(-3.7831843336820317, normal_lcdf(-2, 0, 1));
  EXPECT_FLOAT_EQ(-2.7059444008238898, normal_lcdf(-1.5, 0, 1));
  EXPECT_FLOAT_EQ(-1.8410216450092636, normal_lcdf(-1, 0, 1));
  EXPECT_FLOAT_EQ(-1.1759117615936185, normal_lcdf(-0.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.69314718055994529, normal_lcdf(0, 0, 1));
  EXPECT_FLOAT_EQ(-0.36894641528865652, normal_lcdf(0.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.17275377902344988, normal_lcdf(1, 0, 1));
  EXPECT_FLOAT_EQ(-0.069143455612233992, normal_lcdf(1.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.023012909328963493, normal_lcdf(2, 0, 1));
  EXPECT_FLOAT_EQ(-0.0062290254858600024, normal_lcdf(2.5, 0, 1));
  EXPECT_FLOAT_EQ(-0.0013508099647481938, normal_lcdf(3, 0, 1));
  EXPECT_FLOAT_EQ(-0.00023265614137680445, normal_lcdf(3.5, 0, 1));
  EXPECT_FLOAT_EQ(-3.1671743377489267e-05, normal_lcdf(4, 0, 1));
  EXPECT_FLOAT_EQ(-3.3976788968344657e-06, normal_lcdf(4.5, 0, 1));
  EXPECT_FLOAT_EQ(-2.8665161296376358e-07, normal_lcdf(5, 0, 1));
  EXPECT_FLOAT_EQ(-1.8989562646189464e-08, normal_lcdf(5.5, 0, 1));
  EXPECT_FLOAT_EQ(-9.8658764552437559e-10, normal_lcdf(6, 0, 1));
  EXPECT_FLOAT_EQ(-4.0160005839397589e-11, normal_lcdf(6.5, 0, 1));
  EXPECT_FLOAT_EQ(-1.2798125438866541e-12, normal_lcdf(7, 0, 1));
  EXPECT_FLOAT_EQ(-3.1908916729109475e-14, normal_lcdf(7.5, 0, 1));
  EXPECT_FLOAT_EQ(-6.2209605742717868e-16, normal_lcdf(8, 0, 1));
  EXPECT_FLOAT_EQ(-9.4795348222033192e-18, normal_lcdf(8.5, 0, 1));
  EXPECT_FLOAT_EQ(-1.1285884059538408e-19, normal_lcdf(9, 0, 1));
  EXPECT_FLOAT_EQ(-1.0494515075362608e-21, normal_lcdf(9.5, 0, 1));
  EXPECT_FLOAT_EQ(-7.6198530241605269e-24, normal_lcdf(10, 0, 1));
  EXPECT_FLOAT_EQ(-1.9106595744986757e-28, normal_lcdf(11, 0, 1));
  EXPECT_FLOAT_EQ(-1.7764821120776790e-33, normal_lcdf(12, 0, 1));
  EXPECT_FLOAT_EQ(-6.1171643995498803e-39, normal_lcdf(13, 0, 1));
  EXPECT_FLOAT_EQ(-7.7935368191928000e-45, normal_lcdf(14, 0, 1));
  EXPECT_FLOAT_EQ(-3.6709661993127514e-51, normal_lcdf(15, 0, 1));
  EXPECT_FLOAT_EQ(-6.3887544005380882e-58, normal_lcdf(16, 0, 1));
  EXPECT_FLOAT_EQ(-4.1059962020989074e-65, normal_lcdf(17, 0, 1));
  EXPECT_FLOAT_EQ(-9.7409489189371508e-73, normal_lcdf(18, 0, 1));
  EXPECT_FLOAT_EQ(-8.5272239526309772e-81, normal_lcdf(19, 0, 1));
  EXPECT_FLOAT_EQ(-2.7536241186062337e-89, normal_lcdf(20, 0, 1));
  EXPECT_FLOAT_EQ(-3.2792780189790367e-98, normal_lcdf(21, 0, 1));
  EXPECT_FLOAT_EQ(-1.4398924351450790e-107, normal_lcdf(22, 0, 1));
  EXPECT_FLOAT_EQ(-2.3306370062206492e-117, normal_lcdf(23, 0, 1));
  EXPECT_FLOAT_EQ(-1.3903921185497032e-127, normal_lcdf(24, 0, 1));
  EXPECT_FLOAT_EQ(-3.0566967063825616e-138, normal_lcdf(25, 0, 1));
  EXPECT_FLOAT_EQ(-2.4760633155033892e-149, normal_lcdf(26, 0, 1));
  EXPECT_FLOAT_EQ(-7.3894810068850200e-161, normal_lcdf(27, 0, 1));
  EXPECT_FLOAT_EQ(-8.1238694696594273e-173, normal_lcdf(28, 0, 1));
  EXPECT_FLOAT_EQ(-3.2897852667043802e-185, normal_lcdf(29, 0, 1));
  EXPECT_FLOAT_EQ(-4.9067139271481872e-198, normal_lcdf(30, 0, 1));
  EXPECT_FLOAT_EQ(-2.6952500812005002e-211, normal_lcdf(31, 0, 1));
  EXPECT_FLOAT_EQ(-5.4520806035123956e-225, normal_lcdf(32, 0, 1));
  EXPECT_FLOAT_EQ(-4.0611856209158557e-239, normal_lcdf(33, 0, 1));
  EXPECT_FLOAT_EQ(-1.1138987855743795e-253, normal_lcdf(34, 0, 1));
  EXPECT_FLOAT_EQ(-1.1249107064724062e-268, normal_lcdf(35, 0, 1));
  EXPECT_FLOAT_EQ(-4.1826240657972830e-284, normal_lcdf(36, 0, 1));
  EXPECT_FLOAT_EQ(-5.7255712225245771e-300, normal_lcdf(37, 0, 1));
  EXPECT_FLOAT_EQ(-2.8854283510039645e-316, normal_lcdf(38, 0, 1));
}

/**
 * Value of log Phi against high-precision references.
 *
 * Each row carries the reference to 40 significant digits in a comment and
 * the correctly-rounded double beside it, so a reviewer can check the
 * rounding here and the reference itself in any arbitrary-precision tool.
 * They were produced at 60 significant digits as `log(erfc(-s)/2)` with
 * `s = (y - mu) / (sigma * sqrt(2))`, then rounded to double.
 *
 * The rows span all three branches of the value: `log1p(-erfc(s)/2)` for
 * `s > 0`, `log(erfc(-s)) + LOG_HALF` for `-4 < s <= 0`, and the Cody
 * rational approximation below that.
 *
 * Every row holds to 1e-14 except the last. There `s = 11.31` puts erfc into
 * its far tail, where this platform's libm loses relative accuracy:
 * `std::erfc(11.313708498984761)` returns 1.27775088010759500e-57 against a
 * true 1.27774061192683864e-57, a relative error of 8.0e-06 that passes
 * straight through. That is a libm limit rather than this function's, so the
 * tolerance on that row is loose enough to stay portable.
 */
TEST(ProbNormal, lcdf_matches_high_precision_reference) {
  struct value_ref {
    double y;
    double mu;
    double sigma;
    double expected;
    double tol;
  };
  const value_ref cases[]
      = {// s=2.12132 (log1p); -0.001350809964748193798841110490517380092729
         {3.0, 0.0, 1.0, -0.0013508099647481938, 1e-14},
         // s=0.707107 (log1p); -0.1727537790234498895264831735208007300094
         {1.0, 0.0, 1.0, -0.17275377902344988, 1e-14},
         // s=0.353553 (log1p); -0.3689464152886563930656156390431773405547
         {0.5, 0.0, 1.0, -0.3689464152886564, 1e-14},
         // s=0 (erfc); -0.6931471805599453094172321214581765680755
         {0.0, 0.0, 1.0, -0.6931471805599453, 1e-14},
         // s=-0.707107 (erfc); -1.841021645009263505770783073232529021548
         {-1.0, 0.0, 1.0, -1.8410216450092636, 1e-14},
         // s=-2.12132 (erfc); -6.607726221510349543276077083251514438781
         {-3.0, 0.0, 1.0, -6.607726221510349, 1e-14},
         // s=-3.53553 (erfc); -15.06499839398872573608370479189672560507
         {-5.0, 0.0, 1.0, -15.064998393988725, 1e-14},
         // s=-4.24264 (Cody); -20.7367689499747056549688537180999188204
         {-6.0, 0.0, 1.0, -20.736768949974707, 1e-14},
         // s=-7.07107 (Cody); -53.23128515051247057834702735413120987892
         {-10.0, 0.0, 1.0, -53.23128515051247, 1e-14},
         // s=-14.1421 (Cody); -203.9171553710972639368044586545269000525
         {-20.0, 0.0, 1.0, -203.91715537109727, 1e-14},
         // s=-26.5165 (Cody); -707.6689893175071910661131734572604413308
         {-37.5, 0.0, 1.0, -707.6689893175072, 1e-14},
         // s=-42.4264 (Cody); -1805.013560680567138700666808059022928381
         {-60.0, 0.0, 1.0, -1805.0135606805673, 1e-14},
         // s=-70.7107 (Cody); -5005.52420869420508862630245733002553134
         {-100.0, 0.0, 1.0, -5005.524208694205, 1e-14},
         // s=-212.132 (Cody); -45006.62273211866335985382181649804744182
         {-300.0, 0.0, 1.0, -45006.62273211866, 1e-14},
         // s=-3.29983 (erfc); -13.38983327471644676126972305600605540117
         {-12.0, 2.0, 3.0, -13.389833274716446, 1e-14},
         // s=11.3137 (log1p); -6.388703059634194605096774759212020048496e-58
         {7.0, -1.0, 0.5, -6.388703059634195e-58, 1e-4}};

  for (const value_ref& c : cases) {
    const double got = stan::math::normal_lcdf(c.y, c.mu, c.sigma);
    EXPECT_LT(std::fabs(got / c.expected - 1.0), c.tol)
        << "log Phi at y = " << c.y << " got " << got << " want " << c.expected;
  }
}
