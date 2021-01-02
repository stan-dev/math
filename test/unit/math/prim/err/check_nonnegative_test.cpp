#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>
#include <string>

TEST(ErrorHandlingArr, CheckNonnegativeVectorized) {
  using stan::math::check_nonnegative;
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

TEST(ErrorHandlingArr, CheckNonnegativeVectorized_one_indexed_message) {
  using stan::math::check_nonnegative;
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

TEST(ErrorHandlingArr, CheckNonnegative_nan) {
  using stan::math::check_nonnegative;
  const char* function = "check_nonnegative";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error);
    x[i] = i;
  }
}

TEST(ErrorHandlingScalar, CheckNonnegative) {
  using stan::math::check_nonnegative;
  const char* function = "check_nonnegative";
  double x = 0;

  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with finite x: " << x;

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_nonnegative(function, "x", x))
      << "check_nonnegative should be true with x = Inf: " << x;

  x = -0.01;
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = " << x;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception with x = -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error)
      << "check_nonnegative should throw exception on NaN: " << x;
}

TEST(ErrorHandlingScalar, CheckNonnegative_nan) {
  using stan::math::check_nonnegative;
  const char* function = "check_nonnegative";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_nonnegative(function, "x", nan), std::domain_error);
}
