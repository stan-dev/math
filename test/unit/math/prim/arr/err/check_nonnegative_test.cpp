#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>
#include <string>

using stan::math::check_nonnegative;

TEST(ErrorHandlingScalar, CheckNonnegativeVectorized) {
  int N = 5;
  const char* function = "check_nonnegative";
  std::vector<double> x(N);

  x.assign(N, 0);
  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative(vector) should be true with finite x: " << x[0];

  x.assign(N, std::numeric_limits<double>::infinity());
  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative(vector) should be true with x = Inf: " << x[0];

  x.assign(N, -0.01);
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = " << x[0];

  x.assign(N, -std::numeric_limits<double>::infinity());
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative(vector) should throw an exception with x = -Inf: "
      << x[0];

  x.assign(N, std::numeric_limits<double>::quiet_NaN());
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative(vector) should throw exception on NaN: " << x[0];
}

TEST(ErrorHandlingScalar, CheckNonnegativeVectorized_one_indexed_message) {
  int N = 5;
  const char* function = "check_nonnegative";
  std::vector<double> x(N);
  std::string message;

  x.assign(N, 0);
  x[2] = -1;
  try {
    check_nonnegative(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]"));
}

TEST(ErrorHandlingScalar, CheckNonnegative_nan) {
  const char* function = "check_nonnegative";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error);
    x[i] = i;
  }
}
