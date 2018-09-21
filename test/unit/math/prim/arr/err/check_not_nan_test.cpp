#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <string>

using stan::math::check_not_nan;

TEST(ErrorHandlingScalar, CheckNotNanVectorized) {
  int N = 5;
  const char* function = "check_not_nan";
  std::vector<double> x(N);

  x.assign(N, 0);
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan(vector) should be true with finite x: " << x[0];

  x.assign(N, std::numeric_limits<double>::infinity());
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan(vector) should be true with x = Inf: " << x[0];

  x.assign(N, -std::numeric_limits<double>::infinity());
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan(vector) should be true with x = -Inf: " << x[0];

  x.assign(N, std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(check_not_nan(function, "x", x), std::domain_error)
      << "check_not_nan(vector) should throw exception on NaN: " << x[0];
}

TEST(ErrorHandlingScalar, CheckNotNanVectorized_one_indexed_message) {
  int N = 5;
  const char* function = "check_not_nan";
  std::vector<double> x(N);
  std::string message;

  x.assign(N, 0);
  x[2] = std::numeric_limits<double>::quiet_NaN();
  try {
    check_not_nan(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]")) << message;
}
