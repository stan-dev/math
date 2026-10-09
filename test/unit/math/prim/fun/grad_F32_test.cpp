#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

// converge
TEST(MathPrimScalFun, grad_F32_converges_by_z) {
  std::vector<double> g(6);
  g[0] = 2.290726829685388;
  g[1] = 2.290726829685388;
  g[2] = 2.290726829685388;
  g[3] = -2.290726829685388;
  g[4] = -2.290726829685388;
  g[5] = 6.249999999999999;
  double g_calc[6];
  stan::math::grad_F32(g_calc, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 1e-10);
  for (int i = 0; i < 6; ++i)
    EXPECT_NEAR(g[i], g_calc[i], 1e-8);
}

// FIXME: Can't run till we use a continuation, crosses pole in digamma
// TEST(MathPrimScalFun, grad_F32_polynomial) {
// terminate by zero numerator, no sign-flip
//  EXPECT_NEAR(11.28855722705942,
//              stan::math::grad_F32(1.0, 31.0, -27.0, 19.0, -41.0,
//                .99999, 1e-10), 1e-8);
// }

// FIXME: Can't run till we use a continuation, crosses pole in digamma
// TEST(MathPrimScalFun, grad_F32_short_polynomial) {
// terminate by zero numerator, single-step, no sign-flip
//  std::vector<double> g(6);
//  g[0] = -1.08;
//  g[1] = -0.09;
//  g[2] = -0.07784317107588729;
//  g[3] =  0.108;
//  g[4] =  1.08;
//  g[5] = -1.2;
//  double g_calc[6];
//  stan::math::grad_F32(g_calc, 1.0, 12.0, -1.0, 10.0, 1.0, .9, 1e-20);
//  for (int i=0; i < 6; ++i) {
//    std::cout << i << std::endl;
//    EXPECT_NEAR(g[i], g_calc[i], 1e-8);
//  }
//
// }

// terminate by zero numerator at k = 1; the denominator (b2)_k is zero only
// from k = 2 on, so the polynomial is defined. The derivative with respect
// to a3 is not checked: a3 is the integer that ends the series.
TEST(MathPrimScalFun, grad_F32_short_polynomial_denominator_equal_numerator) {
  std::vector<double> g = {1.2, 0.1, 0.0, -0.12, 1.2, 1.2};
  double g_calc[6];
  stan::math::grad_F32(g_calc, 1.0, 12.0, -1.0, 10.0, -1.0, 1.0);
  for (int i : {0, 1, 3, 4, 5}) {
    EXPECT_NEAR(g[i], g_calc[i], 1e-14);
  }
}

// at pole, should throw
TEST(MathPrimScalFun, grad_F32_short_polynomial_undef) {
  double g_calc[6];
  EXPECT_THROW(stan::math::grad_F32(g_calc, 1.0, 12.0, -2.0, 10.0, -1.0, 1.0),
               std::domain_error);
}

// FIXME: Can't run until we use a continuation, crosses pole in digamma
// TEST(MathPrimScalFun, grad_F32_sign_flip_numerator) {
// converge, single sign flip via numerator
//  std::vector<double> g(6);
//  g[0] = -1.08;
//  g[1] = -0.09;
//  g[2] = -0.07784317107588729;
//  g[3] =  0.108;
//  g[4] =  1.08;
//  g[5] = -1.2;
//  double g_calc[6];
//  stan::math::grad_F32(g_calc, 1.0, -.5, 2.0, 10.0, 1.0, 0.3, 1e-10);
//  for (int i=0; i < 6; ++i)
//    EXPECT_NEAR(g[i], g_calc[i], 1e-8);
//
// }

TEST(MathPrimScalFun, grad_F32_diverge_by_z) {
  // This should throw
  // (Mathematica claims the answer is -10 but it's by continuation
  double g_calc[6];
  EXPECT_THROW(stan::math::grad_F32(g_calc, 1.0, 12.0, 1.0, 10.0, 1.0, 1.1),
               std::domain_error);
}

// convergence, double sign flip
TEST(MathPrimScalFun, grad_F32_double_sign_flip_1) {
  double g_calc[6];
  std::vector<double> g(6);
  g[0] = 0.03692914737912187;
  g[1] = -0.0749983829996620;
  g[2] = -0.0145982785371077;
  g[3] = -0.00367744563081578;
  g[4] = -0.0369291473791218;
  g[5] = 0.1224673892640317;
  stan::math::grad_F32(g_calc, 1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10);
  for (int i = 0; i < 6; ++i)
    EXPECT_NEAR(g[i], g_calc[i], 1e-8);
}

// convergence, double sign flip
TEST(MathPrimScalFun, grad_F32_double_sign_flip_2) {
  double g_calc[6];
  std::vector<double> g(6);
  g[0] = 0.06517384196658483;
  g[1] = -0.1349675922478992;
  g[2] = -0.01422583581177123;
  g[3] = -0.00645591019603955;
  g[4] = -0.0651738419665848;
  g[5] = 0.2147503542109778;
  stan::math::grad_F32(g_calc, 1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10);
  for (int i = 0; i < 6; ++i)
    EXPECT_NEAR(g[i], g_calc[i], 1e-8);
}

// The tests below use z = 1 and terminate by zero numerator at k = 6. The
// reference values are the derivatives of the exact finite sums. The
// derivative with respect to a3 is not checked: a3 is the integer that ends
// the series.

// b2 = a3, so the denominator (b2)_k is zero only from k = 7 on
TEST(MathPrimScalFun, grad_F32_polynomial_denominator_equal_numerator) {
  std::vector<double> g
      = {26.180979740203373,  6.4919047619047619, 0.0,
         -8.0448964404555225, 15.921514834449405, 52.373725818452381};
  double g_calc[6];
  stan::math::grad_F32(g_calc, 1.0, 6.5, -6.0, 5.0, -6.0, 1.0);
  for (int i : {0, 1, 3, 4, 5}) {
    EXPECT_NEAR(g[i], g_calc[i], 1e-12 * std::fabs(g[i]));
  }
}

// a2 and b2 are large
TEST(MathPrimScalFun, grad_F32_polynomial_large_parameters) {
  std::vector<double> g
      = {29.010777844568475, 0.16148712923118207, 0.0,
         -8.114218400390939, 0.32571439023403243, 49.018778179540667};
  double g_calc[6];
  stan::math::grad_F32(g_calc, 1.0, 302.30333970, -6.0, 5.0, -151.75306521,
                       1.0);
  for (int i : {0, 1, 3, 4, 5}) {
    EXPECT_NEAR(g[i], g_calc[i], 1e-12 * std::fabs(g[i]));
  }
}

// a1 = -1 ends the series at k = 1, before (b1)_k is zero from k = 3 on; the
// larger non-positive integer a2 = -5 does not matter. The series is
// 1 + a1 a2 a3 z / (b1 b2). The derivative with respect to a1 is not
// checked: a1 is the integer that ends the series.
TEST(MathPrimScalFun, grad_F32_polynomial_smallest_numerator_ends) {
  std::vector<double> g = {0.0, 0.25, -1.25, -0.625, 1.25, -2.5};
  double g_calc[6];
  stan::math::grad_F32(g_calc, -1.0, -5.0, 1.0, -2.0, 1.0, 0.5);
  for (int i : {1, 2, 3, 4, 5}) {
    EXPECT_NEAR(g[i], g_calc[i], 1e-15);
  }
}

//
// m = {
//  {1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 1e-10},
//  {1.0, 31.0, -27.0, 19.0, -41.0, .99999, 1e-10},
//  {1.0, 12.0, -1.0, 10.0, 1.0, .9},
//  {1.0, 12.0, -1.0, 10.0, -1.0, 1.0},
//  {1.0, -.5, 2.0, 10.0, 1.0, 0.3, 1e-10},
//  {1.0, 12.0, 1.0, 10.0, 1.0, 1.1},
//  {1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10},
//  {1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10}
// }
//
// {2.290726829685388, 2.290726829685388, 2.290726829685388,
// -2.290726829685388, -2.290726829685388, 6.249999999999999},
// IGNORE NO GRADIENT => {22.95673741899624, 1.73991856855972,
// Indeterminate, -2.692768015657851, 1.606385109522632, 59.65366361215116},
// {-1.08, -0.09, -0.07784317107588729, 0.108, 1.08, -1.2},
// THROW => {1.2, 0.1, Indeterminate, -0.12, 1.2, 1.2},
// {-0.03098182751890103, 0.05997687803893869, -0.01554776878589169,
//   0.003126435989341982, 0.03098182751890103, -0.1044310337234595},
// THROW => {-34.02585092994046+31.41592653589793 I,
//            -5.371146427869409+24.22439931425321 I,
//            -34.02585092994046+31.41592653589793 I,
//            7.371146427869409-24.22439931425321 I,
//            34.02585092994046-31.41592653589793 I, 0.},
// {0.03692914737912187, -0.07499838299966203, -0.01459827853710772,
//  -0.003677445630815789, -0.03692914737912187, 0.1224673892640317},
// {0.06517384196658483, -0.1349675922478992, -0.01422583581177123,
//  -0.006455910196039558, -0.06517384196658483, 0.2147503542109778}
//
