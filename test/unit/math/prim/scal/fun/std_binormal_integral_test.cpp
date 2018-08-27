#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/fun/std_binormal_integral.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <unordered_map>
#include <iomanip>

TEST(MathFunctions, binormal_integral_using) {
  using stan::math::std_binormal_integral;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = 0.3;
  double a = nan;
  double b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = 0.3;
  double a = 2;
  double b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = nan;
  double a = 2;
  double b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg) {
  double rho = -1.3;
  double a = 2;
  double b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one) {
  double rho = 1.3;
  double a = 2;
  double b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw) {
  double rho = 0.3;
  double a = 2;
  double b = 1;
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  double rho = 0;
  double a = -0.4;
  double b = 2.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a) * stan::math::Phi(b),
                  stan::math::std_binormal_integral(a, b, rho));

  // Perfectly correlated RVs
  rho = 1;
  a = -3.4;
  b = 3.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a),
                  stan::math::std_binormal_integral(a, b, rho));

  // Perfectly correlated RVs
  rho = 1;
  a = -3.4;
  b = 3.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a),
                  stan::math::std_binormal_integral(a, b, rho));

  // Perfectly anticorrelated RVs
  rho = -1;
  a = 2.4;
  b = 1.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a) + stan::math::Phi(b) - 1,
                  stan::math::std_binormal_integral(a, b, rho));

  // Perfectly anticorrelated RVs
  rho = -1;
  a = -2.4;
  b = 1.7;
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_integral(a, b, rho));

  // a = rho * b
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  EXPECT_FLOAT_EQ(
      0.5 / stan::math::pi() * std::exp(-0.5 * b * b) * std::asin(rho)
          + stan::math::Phi(a) * stan::math::Phi(b),
      stan::math::std_binormal_integral(a, b, rho));

  // b = rho * a
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  EXPECT_FLOAT_EQ(
      0.5 / stan::math::pi() * std::exp(-0.5 * a * a) * std::asin(rho)
          + stan::math::Phi(a) * stan::math::Phi(b),
      stan::math::std_binormal_integral(a, b, rho));
  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1, stan::math::std_binormal_integral(a, b, rho));

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(2.5),
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(0.5),
                  stan::math::std_binormal_integral(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  double rho = 0.3;
  double a = -0.4;
  double b = 2.7;
  EXPECT_FLOAT_EQ(0.344276561500873,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  EXPECT_FLOAT_EQ(0.3445782583896758,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  EXPECT_FLOAT_EQ(0.9937227710497979,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  EXPECT_FLOAT_EQ(0.9997643606337163,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  EXPECT_FLOAT_EQ(2.146032113348184e-06,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = 10;
  EXPECT_FLOAT_EQ(3.397673124738709e-06,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  EXPECT_FLOAT_EQ(4.016000583859118e-11,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = -10;
  EXPECT_NEAR(0, stan::math::std_binormal_integral(a, b, rho), 1e-8);

  rho = 0.99;
  a = -4.5;
  b = -10;
  EXPECT_FLOAT_EQ(7.619853024160583e-24,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.5;
  a = -4.5;
  b = -10;
  EXPECT_FLOAT_EQ(5.612932952882046e-24,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.5;
  a = -4.5;
  b = -4.5;
  EXPECT_FLOAT_EQ(2.348555308541017e-08,
                  stan::math::std_binormal_integral(a, b, rho));
  //  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho),
  //  std::domain_error);

  rho = 0.5;  // integral doesn't converge
  a = -2.5;
  b = -2.5;
  EXPECT_FLOAT_EQ(0.0006693647475263115589547,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.5;
  a = -3.5;
  b = -3.5;
  EXPECT_FLOAT_EQ(8.061334595291052e-14,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.3;
  a = 2.3;
  b = -3.4;
  EXPECT_FLOAT_EQ(std::log(0.0003017720314683548),
                  std::log(stan::math::std_binormal_integral(a, b, rho)));

  rho = 0.3;
  a = -2.5;
  b = -2.0;
  EXPECT_FLOAT_EQ(std::log(0.0007103359768977393L),
                  std::log(stan::math::std_binormal_integral(a, b, rho)));

  rho = 0.7;
  a = -2.5;
  b = -2.0;
  EXPECT_FLOAT_EQ(std::log(0.00301048502825429L),
                  log(stan::math::std_binormal_integral(a, b, rho)));

  rho = 0.9;
  a = -2.5;
  b = -2.0;
  EXPECT_FLOAT_EQ(0.005332931061066137,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.99;
  a = -2.5;
  b = -2.0;
  EXPECT_FLOAT_EQ(0.006209439620767255,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = -2.5;
  b = -2.0;
  EXPECT_FLOAT_EQ(0.006209665325776138,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = -5.5;
  b = -5.0;
  EXPECT_FLOAT_EQ(std::log(1.898956246588779e-8L),
                  std::log(stan::math::std_binormal_integral(a, b, rho)));

  rho = 0.999;
  a = -5.5;
  b = -5.5;
  EXPECT_FLOAT_EQ(std::log(1.707277832561847e-8L),
                  std::log(stan::math::std_binormal_integral(a, b, rho)));

  rho = 0.999;
  a = -15.5;
  b = -15.5;
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.9;  // problem for integrator
  a = -5.5;
  b = -5.5;
  EXPECT_NEAR(1.725633230170963e-30,
              stan::math::std_binormal_integral(a, b, rho), 1e-29);

  rho = -0.9;  // problem for integrator
  a = -7.5;
  b = -7.5;
  EXPECT_NEAR(1.237626803691678e-41,
              stan::math::std_binormal_integral(a, b, rho), 1e-40);

  rho = 0.9;
  a = -5.5;
  b = -5.5;
  EXPECT_FLOAT_EQ(3.675288303261435e-9,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.7;
  a = -0.5;
  b = -0.5;
  EXPECT_FLOAT_EQ(0.0151520415154598,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 1.5;
  b = 1.5;
  EXPECT_FLOAT_EQ(0.866385597462284,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = 1.5;
  b = 1.5;
  EXPECT_FLOAT_EQ(0.930882284854364,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = 5.5;
  b = 5.5;
  EXPECT_FLOAT_EQ(0.999999979093653,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 5.5;
  b = 5.5;
  EXPECT_FLOAT_EQ(0.999999962020875,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 15.5;
  b = 15.5;
  EXPECT_FLOAT_EQ(1.0, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 15.5;
  b = -15.5;
  EXPECT_FLOAT_EQ(0.0, stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = -15.5;
  b = 15.5;
  EXPECT_FLOAT_EQ(0.0, stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = -15.5;
  b = 15.5;
  EXPECT_FLOAT_EQ(1.734460791793871e-54,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = 7.5;
  b = 7.5;
  EXPECT_FLOAT_EQ(0.999999999999964,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 7.5;
  b = 7.5;
  EXPECT_FLOAT_EQ(0.999999999999936,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = -7.5;
  b = 7.5;
  EXPECT_FLOAT_EQ(3.190891672910918e-14,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = -7.5;
  b = 7.5;
  EXPECT_FLOAT_EQ(4.323210067919746e-15,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = 7.5;
  b = -7.5;
  EXPECT_FLOAT_EQ(3.190891672910918e-14,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = 7.5;
  b = -7.5;
  EXPECT_FLOAT_EQ(4.323210067919746e-15,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = 0.999;
  a = -7.5;
  b = -7.5;
  EXPECT_FLOAT_EQ(2.758570666118945e-14,
                  stan::math::std_binormal_integral(a, b, rho));

  rho = -0.999;
  a = -7.5;
  b = -7.5;
  EXPECT_NEAR(0.0, stan::math::std_binormal_integral(a, b, rho), 1e-40);

  // 000
  rho = 0.999;
  a = 9.5;
  b = 9.5;
  EXPECT_FLOAT_EQ(1, stan::math::std_binormal_integral(a, b, rho));

  // 100
  rho = -0.999;
  a = 9.5;
  b = 9.5;
  EXPECT_FLOAT_EQ(1, stan::math::std_binormal_integral(a, b, rho));

  // 010
  rho = 0.999;
  a = -9.5;
  b = 9.5;
  EXPECT_FLOAT_EQ(1.049451507536276e-21,
                  stan::math::std_binormal_integral(a, b, rho));

  // 110
  rho = -0.999;
  a = -9.5;
  b = 9.5;
  EXPECT_FLOAT_EQ(1.784740987873031e-22,
                  stan::math::std_binormal_integral(a, b, rho));

  // 011
  rho = 0.999;
  a = -9.5;
  b = -9.5;
  EXPECT_FLOAT_EQ(8.70977408748972e-22,
                  stan::math::std_binormal_integral(a, b, rho));

  // 111
  rho = -0.999;
  a = -9.5;
  b = -9.5;
  EXPECT_FLOAT_EQ(0.0, stan::math::std_binormal_integral(a, b, rho));

  // 010
  rho = 0.999;
  a = -35.5;
  b = 35.5;
  EXPECT_FLOAT_EQ(-634.6142631550883379532,
                  std::log(stan::math::std_binormal_integral(a, b, rho)));

  // 110
  rho = -0.999;
  a = -35.5;
  b = 35.5;
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_integral(a, b, rho));

  // 011
  rho = 0.999;
  a = -35.5;
  b = -35.5;
  EXPECT_FLOAT_EQ(2.457691540662217e-276,
                  stan::math::std_binormal_integral(a, b, rho));

  // 111
  rho = -0.999;
  a = -35.5;
  b = -35.5;
  EXPECT_FLOAT_EQ(0.0, stan::math::std_binormal_integral(a, b, rho));

  rho = 0.9999;
  a = -25;
  b = -25;
  for (int i = 0; i < 198; ++i) {
    for (int j = 0; j < 198; ++j) {
      for (int k = 0; k < 198; ++k) {
        EXPECT_GE(stan::math::std_binormal_integral(a, b, rho), 0)
            << "a: " << a << " b: " << b << " r: " << rho << std::endl;
        b += 2.0 * 25.0 / 197.0;
      }
      a += 2.0 * 25.0 / 197.0;
    }
    rho -= 0.01;
  }
}
