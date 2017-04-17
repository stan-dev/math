#include <stan/math/prim/scal.hpp>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, atanh) {
  using stan::math::atanh;
  EXPECT_FLOAT_EQ(-3.341680472883148, atanh(-0.9975));
  EXPECT_FLOAT_EQ(0, atanh(0));
  EXPECT_FLOAT_EQ(0.5493061443340549, atanh(0.5));
}

TEST(MathFunctions, atanh_exception) {
  using stan::math::atanh;
  EXPECT_THROW(atanh(-2), std::domain_error);
  EXPECT_THROW(atanh(2.0), std::domain_error);
}

TEST(MathFunctions, atanh_inf_return) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            stan::math::atanh(-1));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            stan::math::atanh(1));
}

TEST(MathFunctions, atanh_nan) {
  using stan::math::atanh;
  EXPECT_PRED1(boost::math::isnan<double>,
               stan::math::atanh(std::numeric_limits<double>::quiet_NaN()));
}
