#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/util.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

TEST(AgradRev,if_else) {
  using stan::math::var;
  using stan::math::if_else;
  using stan::math::if_else;
  
  EXPECT_FLOAT_EQ(1.0,if_else(true,var(1.0),var(2.0)).val());
  EXPECT_FLOAT_EQ(2.0,if_else(false,var(1.0),var(2.0)).val());

  EXPECT_FLOAT_EQ(1.0,if_else(true,1.0,var(2.0)).val());
  EXPECT_FLOAT_EQ(2.0,if_else(false,1.0,var(2.0)).val());

  EXPECT_FLOAT_EQ(1.0,if_else(true,var(1.0),2.0).val());
  EXPECT_FLOAT_EQ(2.0,if_else(false,var(1.0),2.0).val());
}

TEST(AgradRev, if_else_nan) {
  using stan::math::if_else;

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::math::var nan_v = std::numeric_limits<double>::quiet_NaN();
  stan::math::var a_v = 1.2;

  EXPECT_FLOAT_EQ(1.2, if_else(true, 1.2, nan_v).val());
  EXPECT_FLOAT_EQ(1.2, if_else(true, a_v, nan).val());
  EXPECT_FLOAT_EQ(1.2, if_else(true, a_v, nan_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, 1.2, nan_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, a_v, nan).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, a_v, nan_v).val());

  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan_v, 2.4).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan, a_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan_v, a_v).val());

  a_v = 2.4;
  EXPECT_FLOAT_EQ(2.4, if_else(false, nan_v, 2.4).val());
  EXPECT_FLOAT_EQ(2.4, if_else(false, nan, a_v).val());
  EXPECT_FLOAT_EQ(2.4, if_else(false, nan_v, a_v).val());

  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan_v, nan).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan, nan_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(true, nan_v, nan_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, nan, nan_v).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, nan_v, nan).val());
  EXPECT_PRED1(boost::math::isnan<double>,
               if_else(false, nan_v, nan_v).val());
}

TEST(AgradRev, check_varis_on_stack) {
  stan::math::var x = 1.0;
  stan::math::var y = 2.0;
  test::check_varis_on_stack(stan::math::if_else(true, x, y));
  test::check_varis_on_stack(stan::math::if_else(false, x, y));
  test::check_varis_on_stack(stan::math::if_else(true, x, 2.0));
  test::check_varis_on_stack(stan::math::if_else(false, x, 2.0));
  test::check_varis_on_stack(stan::math::if_else(true, 1.0, y));
  test::check_varis_on_stack(stan::math::if_else(false, 1.0, y));
}
