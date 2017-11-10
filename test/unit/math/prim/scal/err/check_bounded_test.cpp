#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

using stan::math::check_bounded;

TEST(ErrorHandlingScalar,CheckBounded_x) {
  const char* function = "check_bounded";
  const char* name = "x";
  double x = 0;
  double low = -1;
  double high = 1;
 
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high))
    << "check_bounded should be TRUE with x: " << x << " and bounds: " << low << ", " << high;
  
  x = low;
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high)) 
    << "check_bounded should be TRUE with x: " << x << " equal to the lower bound: " << low;

  x = high;
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high)) 
    << "check_bounded should be TRUE with x: " << x << " equal to the lower bound: " << low;

  x = low-1;
  EXPECT_THROW(check_bounded(function, name, x, low, high), 
               std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << high << ", " << low;  
  
  x = high+1;
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << high << ", " << low;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << high << ", " << low;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << high << ", " << low;

  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << high << ", " << low;  
}
TEST(ErrorHandlingScalar,CheckBounded_Low) {
  const char* function = "check_bounded";
  const char* name = "x";
  double x = 0;
  double low = -1;
  double high = 1;
 
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high))
    << "check_bounded should be true x: " << x << " and bounds: " << low << ", " << high;
  
  low = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high))
    << "check_bounded should be TRUE with x: " << x << " and bounds: " << low << ", " << high;

  low = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_bounded(function, name, x, low, high), 
               std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << low << ", " << high;

  low = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << low << ", " << high;
}
TEST(ErrorHandlingScalar,CheckBounded_High) {
  const char* function = "check_bounded";
  const char* name = "x";
  double x = 0;
  double low = -1;
  double high = 1;
 
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high))
    << "check_bounded should be true x: " << x << " and bounds: " << low << ", " << high;

  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_bounded(function, name, x, low, high)) 
    << "check_bounded should be TRUE with x: " << x << " and bounds: " << low << ", " << high;
  
  high = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << low << ", " << high;

  high = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_bounded(function, name, x, low, high), std::domain_error) 
    << "check_bounded should throw with x: " << x << " and bounds: " << low << ", " << high;
}
TEST(ErrorHandlingScalar,CheckBounded_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  const char* function = "check_bounded";
  const char* name = "x";
  double x = 0;
  double low = -1;
  double high = 1;

  EXPECT_THROW(check_bounded(function, name, nan, low, high),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, x, nan, high),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, x, low, nan),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, nan, nan, high),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, nan, low, nan),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, x, nan, nan),
               std::domain_error);
  EXPECT_THROW(check_bounded(function, name, nan, nan, nan),
               std::domain_error);
}
