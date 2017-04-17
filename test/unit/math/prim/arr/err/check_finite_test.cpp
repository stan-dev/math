#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

using stan::math::check_finite;


// ---------- check_finite: vector tests ----------
TEST(ErrorHandlingScalar,CheckFinite_Vector) {
  const char* function = "check_finite";
  std::vector<double> x;
  
  x.clear();
  x.push_back (-1);
  x.push_back (0);
  x.push_back (1);
  ASSERT_NO_THROW(check_finite(function, "x", x)) 
    << "check_finite should be true with finite x";

  x.clear();
  x.push_back(-1);
  x.push_back(0);
  x.push_back(std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error) 
    << "check_finite should throw exception on Inf";

  x.clear();
  x.push_back(-1);
  x.push_back(0);
  x.push_back(-std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
    << "check_finite should throw exception on -Inf";
  
  x.clear();
  x.push_back(-1);
  x.push_back(0);
  x.push_back(std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
    << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingScalar,CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x;
  x.push_back (nan);
  x.push_back (0);
  x.push_back (1);
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x[0] = 1.0;
  x[1] = nan;
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x[1] = 0.0;
  x[2] = nan;
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);
}
