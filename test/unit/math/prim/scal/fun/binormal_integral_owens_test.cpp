#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/fun/binormal_integral_owens.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <unordered_map>

TEST(MathFunctions, binormal_integral_using) { using stan::math::binormal_integral_owens; }

TEST(MathFunctions, binormal_integral_throw_RV_1_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = 0.3;
  double a = nan;
  double b = 2;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = 0.3;
  double a = 2;
  double b = nan;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  double rho = nan;
  double a = 2;
  double b = 2;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg) {
  double rho = -1.3;
  double a = 2;
  double b = 1;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one) {
  double rho = 1.3;
  double a = 2;
  double b = 1;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw) {
  double rho = 0.3;
  double a = 2;
  double b = 1;
  EXPECT_NO_THROW(stan::math::binormal_integral_owens(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  double rho = 0;
  double a = -0.4;
  double b = 2.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a) * stan::math::Phi(b), stan::math::binormal_integral_owens(a, b, rho));

  // Perfectly correlated RVs 
  rho = 1;
  a = -3.4;
  b = 3.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a), stan::math::binormal_integral_owens(a, b, rho));
  
  // Perfectly correlated RVs 
  rho = 1;
  a = -3.4;
  b = 3.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a), stan::math::binormal_integral_owens(a, b, rho));

  // Perfectly anticorrelated RVs 
  rho = -1;
  a = 2.4;
  b = 1.7;
  EXPECT_FLOAT_EQ(stan::math::Phi(a) + stan::math::Phi(b) - 1, stan::math::binormal_integral_owens(a, b, rho));
  
  // Perfectly anticorrelated RVs 
  rho = -1;
  a = -2.4;
  b = 1.7;
  EXPECT_FLOAT_EQ(0, stan::math::binormal_integral_owens(a, b, rho));
  
  // a = rho * b 
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi()
                  * std::exp(-0.5 * b * b)
                  * std::asin(rho)
                  + stan::math::Phi(a) * stan::math::Phi(b), 
                  stan::math::binormal_integral_owens(a, b, rho));
  
  // b = rho * a 
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi()
                  * std::exp(-0.5 * a * a)
                  * std::asin(rho)
                  + stan::math::Phi(a) * stan::math::Phi(b), 
                  stan::math::binormal_integral_owens(a, b, rho));
  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1,stan::math::binormal_integral_owens(a, b, rho)); 

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(0,stan::math::binormal_integral_owens(a, b, rho)); 

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(0,stan::math::binormal_integral_owens(a, b, rho)); 

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(1,stan::math::binormal_integral_owens(a, b, rho)); 

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),stan::math::binormal_integral_owens(a, b, rho)); 

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),stan::math::binormal_integral_owens(a, b, rho)); 

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(2.5),stan::math::binormal_integral_owens(a, b, rho)); 

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  EXPECT_FLOAT_EQ(stan::math::Phi(0.5),stan::math::binormal_integral_owens(a, b, rho)); 
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  double rho = 0.3;
  double a = -0.4;
  double b = 2.7;
  EXPECT_FLOAT_EQ(0.344276561500873, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  EXPECT_FLOAT_EQ(0.3445782583896758, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  EXPECT_FLOAT_EQ(0.9937227710497979, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  EXPECT_FLOAT_EQ(0.9997643606337163, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  EXPECT_FLOAT_EQ(2.146032113348184e-06, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = 10;
  EXPECT_FLOAT_EQ(3.397673124738709e-06, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  EXPECT_FLOAT_EQ(4.016000583859118e-11, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.99;
  a = -4.5;
  b = -10;
  EXPECT_FLOAT_EQ(0, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.99;
  a = -4.5;
  b = -10;
  EXPECT_FLOAT_EQ(7.619853024160583e-24, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.5;
  a = -4.5;
  b = -10;
  EXPECT_FLOAT_EQ(5.612932952882069e-24, stan::math::binormal_integral_owens(a, b, rho));

  rho = 0.5;
  a = -4.5;
  b = -4.5;
  EXPECT_FLOAT_EQ(2.348555308541017e-08, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  EXPECT_FLOAT_EQ(2.437354593037188e-08, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  EXPECT_FLOAT_EQ(2.437354593037188e-08, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.5;
  a = -3.5;
  b = -3.5;
  EXPECT_FLOAT_EQ(8.061334595291052e-14, stan::math::binormal_integral_owens(a, b, rho));

  rho = -0.3;
  a = 2.3;
  b = -3.4;
  EXPECT_FLOAT_EQ(std::log(0.0003017720314683548), std::log(stan::math::binormal_integral_owens(a, b, rho)));

  rho = 0.9999;
  a = -15;
  b = -15;
  for (int i = 0; i < 198; ++i) {
    for (int j = 0; j < 198; ++j) {
      for (int k = 0; k < 198; ++k) {
        EXPECT_GE(stan::math::binormal_integral_owens(a, b, rho), 0) << "a: " << a << " b: " << b << " r: " << rho << std::endl;
        b += 14.0 / 197;
      }
      a += 14.0 / 197;
    }
    rho -= 0.01;
  }
}
