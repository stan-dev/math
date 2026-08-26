#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbStdNormal, cdf_log_matches_lcdf) {
  double y = 0.8;

  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf(y)),
                  (stan::math::std_normal_cdf_log(y)));
  EXPECT_FLOAT_EQ((stan::math::std_normal_lcdf<double>(y)),
                  (stan::math::std_normal_cdf_log<double>(y)));
}

TEST(ProbStdNormal, lcdf_tails) {
  using stan::math::std_normal_lcdf;
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
  //       "EXPECT_FLOAT_EQ(%#.17g, std_normal_lcdf(%.17g));\n",
  //       pnorm(q[i], lower.tail = TRUE, log.p = TRUE),
  //       q[i]
  //     )
  //   )
  // }

  EXPECT_FLOAT_EQ(-5000000000000019.0, std_normal_lcdf(-100000000));
  EXPECT_FLOAT_EQ(-50000000000017.039, std_normal_lcdf(-10000000));
  EXPECT_FLOAT_EQ(-500000000014.73444, std_normal_lcdf(-1000000));
  EXPECT_FLOAT_EQ(-5000000012.4318638, std_normal_lcdf(-100000));
  EXPECT_FLOAT_EQ(-50000010.129278913, std_normal_lcdf(-10000));
  EXPECT_FLOAT_EQ(-500007.82669481216, std_normal_lcdf(-1000));
  EXPECT_FLOAT_EQ(-5005.5242086942053, std_normal_lcdf(-100));
  EXPECT_FLOAT_EQ(-1254.8313611394199, std_normal_lcdf(-50));
  EXPECT_FLOAT_EQ(-1205.3111748916654, std_normal_lcdf(-49));
  EXPECT_FLOAT_EQ(-1156.7905731019453, std_normal_lcdf(-48));
  EXPECT_FLOAT_EQ(-1109.2695383172531, std_normal_lcdf(-47));
  EXPECT_FLOAT_EQ(-1062.7480519624305, std_normal_lcdf(-46));
  EXPECT_FLOAT_EQ(-1017.2260942419524, std_normal_lcdf(-45));
  EXPECT_FLOAT_EQ(-972.70364403073665, std_normal_lcdf(-44));
  EXPECT_FLOAT_EQ(-929.18067875247391, std_normal_lcdf(-43));
  EXPECT_FLOAT_EQ(-886.65717424372951, std_normal_lcdf(-42));
  EXPECT_FLOAT_EQ(-845.13310460177456, std_normal_lcdf(-41));
  EXPECT_FLOAT_EQ(-804.60844201375380, std_normal_lcdf(-40));
  EXPECT_FLOAT_EQ(-765.08315656437753, std_normal_lcdf(-39));
  EXPECT_FLOAT_EQ(-726.55721601882010, std_normal_lcdf(-38));
  EXPECT_FLOAT_EQ(-689.03058557689064, std_normal_lcdf(-37));
  EXPECT_FLOAT_EQ(-652.50322759379844, std_normal_lcdf(-36));
  EXPECT_FLOAT_EQ(-616.97510126192253, std_normal_lcdf(-35));
  EXPECT_FLOAT_EQ(-582.44616224687172, std_normal_lcdf(-34));
  EXPECT_FLOAT_EQ(-548.91636226973810, std_normal_lcdf(-33));
  EXPECT_FLOAT_EQ(-516.38564862572537, std_normal_lcdf(-32));
  EXPECT_FLOAT_EQ(-484.85396362717927, std_normal_lcdf(-31));
  EXPECT_FLOAT_EQ(-454.32124395634321, std_normal_lcdf(-30));
  EXPECT_FLOAT_EQ(-424.78741990973015, std_normal_lcdf(-29));
  EXPECT_FLOAT_EQ(-396.25241451163106, std_normal_lcdf(-28));
  EXPECT_FLOAT_EQ(-368.71614246865636, std_normal_lcdf(-27));
  EXPECT_FLOAT_EQ(-342.17850892992783, std_normal_lcdf(-26));
  EXPECT_FLOAT_EQ(-316.63940800802027, std_normal_lcdf(-25));
  EXPECT_FLOAT_EQ(-292.09872100320780, std_normal_lcdf(-24));
  EXPECT_FLOAT_EQ(-268.55631425686312, std_normal_lcdf(-23));
  EXPECT_FLOAT_EQ(-246.01203653738091, std_normal_lcdf(-22));
  EXPECT_FLOAT_EQ(-224.46571583141449, std_normal_lcdf(-21));
  EXPECT_FLOAT_EQ(-203.91715537109727, std_normal_lcdf(-20));
  EXPECT_FLOAT_EQ(-184.36612866916096, std_normal_lcdf(-19));
  EXPECT_FLOAT_EQ(-165.81237325071419, std_normal_lcdf(-18));
  EXPECT_FLOAT_EQ(-148.25558265098039, std_normal_lcdf(-17));
  EXPECT_FLOAT_EQ(-131.69539607375970, std_normal_lcdf(-16));
  EXPECT_FLOAT_EQ(-116.13138484571169, std_normal_lcdf(-15));
  EXPECT_FLOAT_EQ(-101.56303440744996, std_normal_lcdf(-14));
  EXPECT_FLOAT_EQ(-87.989719971022524, std_normal_lcdf(-13));
  EXPECT_FLOAT_EQ(-75.410673001568796, std_normal_lcdf(-12));
  EXPECT_FLOAT_EQ(-63.824934094423718, std_normal_lcdf(-11));
  EXPECT_FLOAT_EQ(-53.231285150512470, std_normal_lcdf(-10));
  EXPECT_FLOAT_EQ(-48.306019298965232, std_normal_lcdf(-9.5));
  EXPECT_FLOAT_EQ(-43.628149113332114, std_normal_lcdf(-9));
  EXPECT_FLOAT_EQ(-39.197396428217672, std_normal_lcdf(-8.5));
  EXPECT_FLOAT_EQ(-35.013437159914552, std_normal_lcdf(-8));
  EXPECT_FLOAT_EQ(-31.075890902890002, std_normal_lcdf(-7.5));
  EXPECT_FLOAT_EQ(-27.384307498811076, std_normal_lcdf(-7));
  EXPECT_FLOAT_EQ(-23.938149495161838, std_normal_lcdf(-6.5));
  EXPECT_FLOAT_EQ(-20.736768949974707, std_normal_lcdf(-6));
  EXPECT_FLOAT_EQ(-17.779376352625260, std_normal_lcdf(-5.5));
  EXPECT_FLOAT_EQ(-15.064998393988725, std_normal_lcdf(-5));
  EXPECT_FLOAT_EQ(-12.592419735713079, std_normal_lcdf(-4.5));
  EXPECT_FLOAT_EQ(-10.360101486527292, std_normal_lcdf(-4));
  EXPECT_FLOAT_EQ(-8.3660653083440941, std_normal_lcdf(-3.5));
  EXPECT_FLOAT_EQ(-6.6077262215103492, std_normal_lcdf(-3));
  EXPECT_FLOAT_EQ(-5.0816482772786902, std_normal_lcdf(-2.5));
  EXPECT_FLOAT_EQ(-3.7831843336820317, std_normal_lcdf(-2));
  EXPECT_FLOAT_EQ(-2.7059444008238898, std_normal_lcdf(-1.5));
  EXPECT_FLOAT_EQ(-1.8410216450092636, std_normal_lcdf(-1));
  EXPECT_FLOAT_EQ(-1.1759117615936185, std_normal_lcdf(-0.5));
  EXPECT_FLOAT_EQ(-0.69314718055994529, std_normal_lcdf(0));
  EXPECT_FLOAT_EQ(-0.36894641528865652, std_normal_lcdf(0.5));
  EXPECT_FLOAT_EQ(-0.17275377902344988, std_normal_lcdf(1));
  EXPECT_FLOAT_EQ(-0.069143455612233992, std_normal_lcdf(1.5));
  EXPECT_FLOAT_EQ(-0.023012909328963493, std_normal_lcdf(2));
  EXPECT_FLOAT_EQ(-0.0062290254858600024, std_normal_lcdf(2.5));
  EXPECT_FLOAT_EQ(-0.0013508099647481938, std_normal_lcdf(3));
  EXPECT_FLOAT_EQ(-0.00023265614137680445, std_normal_lcdf(3.5));
  EXPECT_FLOAT_EQ(-3.1671743377489267e-05, std_normal_lcdf(4));
  EXPECT_FLOAT_EQ(-3.3976788968344657e-06, std_normal_lcdf(4.5));
  EXPECT_FLOAT_EQ(-2.8665161296376358e-07, std_normal_lcdf(5));
  EXPECT_FLOAT_EQ(-1.8989562646189464e-08, std_normal_lcdf(5.5));
  EXPECT_FLOAT_EQ(-9.8658764552437559e-10, std_normal_lcdf(6));
  EXPECT_FLOAT_EQ(-4.0160005839397589e-11, std_normal_lcdf(6.5));
  EXPECT_FLOAT_EQ(-1.2798125438866541e-12, std_normal_lcdf(7));
  EXPECT_FLOAT_EQ(-3.1908916729109475e-14, std_normal_lcdf(7.5));
  EXPECT_FLOAT_EQ(-6.2209605742717868e-16, std_normal_lcdf(8));
  EXPECT_FLOAT_EQ(-9.4795348222033192e-18, std_normal_lcdf(8.5));
  EXPECT_FLOAT_EQ(-1.1285884059538408e-19, std_normal_lcdf(9));
  EXPECT_FLOAT_EQ(-1.0494515075362608e-21, std_normal_lcdf(9.5));
  EXPECT_FLOAT_EQ(-7.6198530241605269e-24, std_normal_lcdf(10));
  EXPECT_FLOAT_EQ(-1.9106595744986757e-28, std_normal_lcdf(11));
  EXPECT_FLOAT_EQ(-1.7764821120776790e-33, std_normal_lcdf(12));
  EXPECT_FLOAT_EQ(-6.1171643995498803e-39, std_normal_lcdf(13));
  EXPECT_FLOAT_EQ(-7.7935368191928000e-45, std_normal_lcdf(14));
  EXPECT_FLOAT_EQ(-3.6709661993127514e-51, std_normal_lcdf(15));
  EXPECT_FLOAT_EQ(-6.3887544005380882e-58, std_normal_lcdf(16));
  EXPECT_FLOAT_EQ(-4.1059962020989074e-65, std_normal_lcdf(17));
  EXPECT_FLOAT_EQ(-9.7409489189371508e-73, std_normal_lcdf(18));
  EXPECT_FLOAT_EQ(-8.5272239526309772e-81, std_normal_lcdf(19));
  EXPECT_FLOAT_EQ(-2.7536241186062337e-89, std_normal_lcdf(20));
  EXPECT_FLOAT_EQ(-3.2792780189790367e-98, std_normal_lcdf(21));
  EXPECT_FLOAT_EQ(-1.4398924351450790e-107, std_normal_lcdf(22));
  EXPECT_FLOAT_EQ(-2.3306370062206492e-117, std_normal_lcdf(23));
  EXPECT_FLOAT_EQ(-1.3903921185497032e-127, std_normal_lcdf(24));
  EXPECT_FLOAT_EQ(-3.0566967063825616e-138, std_normal_lcdf(25));
  EXPECT_FLOAT_EQ(-2.4760633155033892e-149, std_normal_lcdf(26));
  EXPECT_FLOAT_EQ(-7.3894810068850200e-161, std_normal_lcdf(27));
  EXPECT_FLOAT_EQ(-8.1238694696594273e-173, std_normal_lcdf(28));
  EXPECT_FLOAT_EQ(-3.2897852667043802e-185, std_normal_lcdf(29));
  EXPECT_FLOAT_EQ(-4.9067139271481872e-198, std_normal_lcdf(30));
  EXPECT_FLOAT_EQ(-2.6952500812005002e-211, std_normal_lcdf(31));
  EXPECT_FLOAT_EQ(-5.4520806035123956e-225, std_normal_lcdf(32));
  EXPECT_FLOAT_EQ(-4.0611856209158557e-239, std_normal_lcdf(33));
  EXPECT_FLOAT_EQ(-1.1138987855743795e-253, std_normal_lcdf(34));
  EXPECT_FLOAT_EQ(-1.1249107064724062e-268, std_normal_lcdf(35));
  EXPECT_FLOAT_EQ(-4.1826240657972830e-284, std_normal_lcdf(36));
  EXPECT_FLOAT_EQ(-5.7255712225245771e-300, std_normal_lcdf(37));
  EXPECT_FLOAT_EQ(-2.8854283510039645e-316, std_normal_lcdf(38));
}
