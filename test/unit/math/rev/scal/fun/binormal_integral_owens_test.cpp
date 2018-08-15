#include <stan/math/rev/scal.hpp>
#include <stan/math/prim/scal/fun/binormal_integral_owens.hpp>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, binormal_integral_using) { using stan::math::binormal_integral_owens; }

TEST(MathFunctions, binormal_integral_throw_RV_1_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<double>::quiet_NaN();
  var rho = 0.3;
  var a = nan;
  var b = 2;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<var>::quiet_NaN();
  var rho = 0.3;
  var a = 2;
  var b = nan;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_vvv) {
  using stan::math::var;
  var nan = std::numeric_limits<double>::quiet_NaN();
  var rho = nan;
  var a = 2;
  var b = 2;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_vvv) {
  using stan::math::var;
  var rho = -1.3;
  var a = 2;
  var b = 1;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_vvv) {
  using stan::math::var;
  var rho = 1.3;
  var a = 2;
  var b = 1;
  EXPECT_THROW(stan::math::binormal_integral_owens(a, b, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_vvv) {
  using stan::math::var;
  var rho = 0.3;
  var a = 2;
  var b = 1;
  EXPECT_NO_THROW(stan::math::binormal_integral_owens(a, b, rho));
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  using stan::math::var;
  var rho = 0;
  var a = -0.4;
  var b = 2.7;
  var f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()) * stan::math::Phi(b.val()), f.val());

  // Perfectly correlated RVs 
  rho = 1;
  a = -3.4;
  b = 3.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()), f.val());
  
  // Perfectly anticorrelated RVs 
  rho = -1;
  a = 2.4;
  b = 1.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(a.val()) + stan::math::Phi(b.val()) - 1, f.val());
  
  // Perfectly anticorrelated RVs 
  rho = -1;
  a = -2.4;
  b = 1.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0, f.val());
  
  // a = rho * b 
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi()
                  * std::exp(-0.5 * b.val() * b.val())
                  * std::asin(rho.val())
                  + stan::math::Phi(a.val()) * stan::math::Phi(b.val()), 
                  f.val());
  // b = rho * a 
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.5 / stan::math::pi()
                  * std::exp(-0.5 * a.val() * a.val())
                  * std::asin(rho.val())
                  + stan::math::Phi(a.val()) * stan::math::Phi(b.val()), 
                  f.val());
  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(1,f.val()); 

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0,f.val()); 

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0,f.val()); 

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(1,f.val()); 

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),f.val()); 

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(1.5),f.val()); 

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(2.5),f.val()); 

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(stan::math::Phi(0.5),f.val()); 
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  using stan::math::var;
  var rho = 0.3;
  var a = -0.4;
  var b = 2.7;
  var f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.344276561500873, f.val());

  rho = 0.99;
  a = -0.4;
  b = 2.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.3445782583896758, f.val());

  rho = 0.99;
  a = 2.5;
  b = 2.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.9937227710497979, f.val());

  rho = 0.99;
  a = 3.5;
  b = 3.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0.9997643606337163, f.val());

  rho = -0.99;
  a = -4.5;
  b = 4.7;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(2.146032113348184e-06, f.val());

  rho = -0.99;
  a = -4.5;
  b = 10;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(3.397673124738709e-06, f.val());

  rho = 0.99;
  a = 4.5;
  b = -6.5;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(4.016000583859118e-11, f.val());

  rho = -0.99;
  a = -4.5;
  b = -10;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(0, f.val());

  rho = 0.99;
  a = -4.5;
  b = -10;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(7.619853024160583e-24, f.val());

  rho = 0.5;
  a = -4.5;
  b = -10;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(5.612932952882069e-24, f.val());

  rho = 0.5;
  a = -4.5;
  b = -4.5;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(2.348555308541017e-08, f.val());

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val());

  rho = -0.5;
  a = -2.5;
  b = -2.5;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(2.437354593037188e-08, f.val());

  rho = -0.5;
  a = -3.5;
  b = -3.5;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(8.061334595291052e-14, f.val());

  rho = -0.3;
  a = -3.3;
  b = -3.4;
  f = stan::math::binormal_integral_owens(a, b, rho);
  EXPECT_FLOAT_EQ(-21.05454315638917, log(f.val()));
}
