#include <stan/math/prim/scal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, log_modified_bessel_first_kind) {
  using stan::math::log_modified_bessel_first_kind;
  using std::sqrt;
  using std::log;

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(sqrt(3), sqrt(2)),
                  log(boost::math::cyl_bessel_i( sqrt(3), sqrt(2))));

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(0.5, 10),
                  log(sqrt(2 / (M_PI * 10)) * std::sinh(10)));

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(-0.5, 10),
                  log(sqrt(2 / (M_PI * 10)) * std::cosh(10)));

  // code branches for z > 100 and v sufficiently small
  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(.1, 100),
                  log_modified_bessel_first_kind(.1, 100 + 1e-16));

  // code branches when v is -1
  EXPECT_FLOAT_EQ(-0.5706479874908312,
                  log_modified_bessel_first_kind(-1, 1));
  EXPECT_FLOAT_EQ(-0.5706479874908312,
                  log_modified_bessel_first_kind(-1 + 1e-16, 1));

  // integers inputs get promoted correctly
  EXPECT_FLOAT_EQ(0.235914358507179,
                  log_modified_bessel_first_kind(0, 1));

  // limiting cases
  EXPECT_EQ(0, log_modified_bessel_first_kind(0,0));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            log_modified_bessel_first_kind(1,0));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            log_modified_bessel_first_kind(-0.1,0));
}

TEST(MathFunctions, log_modified_bessel_first_kind_throw) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  using stan::math::log_modified_bessel_first_kind;

  EXPECT_THROW(log_modified_bessel_first_kind(1, nan), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(nan, 1), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(.5, -2), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(-2,  1), std::domain_error);
}
