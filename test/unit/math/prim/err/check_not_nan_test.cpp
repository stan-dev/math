#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <string>
#include <stdexcept>

TEST(ErrorHandlingArr, CheckNotNanVectorized) {
  using stan::math::check_not_nan;
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

TEST(ErrorHandlingArr, CheckNotNanVectorized_one_indexed_message) {
  using stan::math::check_not_nan;
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

TEST(ErrorHandlingMatrix, checkNotNanEigenRow) {
  using stan::math::check_not_nan;
  stan::math::vector_d y;
  y.resize(3);
  y << 1, 2, 3;

  EXPECT_NO_THROW(stan::math::check_not_nan("checkNotNanEigenRow(%1)", "y", y));
  EXPECT_NO_THROW(stan::math::check_not_nan("checkNotNanEigenRow(%1)", "y", y));

  y(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::check_not_nan("checkNotNanEigenRow", "y", y),
               std::domain_error);
  EXPECT_THROW(stan::math::check_not_nan("checkNotNanEigenRow", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingScalar, CheckNotNan) {
  using stan::math::check_not_nan;
  const char* function = "check_not_nan";
  double x = 0;

  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan should be true with finite x: " << x;

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan should be true with x = Inf: " << x;

  x = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_not_nan(function, "x", x))
      << "check_not_nan should be true with x = -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_not_nan(function, "x", x), std::domain_error)
      << "check_not_nan should throw exception on NaN: " << x;
}
