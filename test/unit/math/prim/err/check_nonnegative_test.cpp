#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>
#include <string>

using stan::math::check_nonnegative;
const char* function = "check_nonnegative";

TEST(ErrorHandlingArr, CheckNonnegativeVectorized) {
  int N = 5;
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
  int N = 5;
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
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {1, 2, 3};

  for (size_t i = 0; i < x.size(); i++) {
    x[i] = nan;
    EXPECT_THROW(check_nonnegative(function, "x", x), std::domain_error);
    x[i] = i;
  }
}

TEST(ErrorHandlingScalar, CheckNonnegative) {
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
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_nonnegative(function, "x", nan), std::domain_error);
}

TEST(ErrorHandlingScalar, CheckNonnegative_0) {
  EXPECT_NO_THROW(check_nonnegative(function, "x", 0U));
  EXPECT_NO_THROW(check_nonnegative(function, "x", (size_t)0));
  EXPECT_NO_THROW(check_nonnegative(function, "x", 0.0));
  EXPECT_NO_THROW(check_nonnegative(function, "x", -0.0));
  EXPECT_NO_THROW(check_nonnegative(function, "x", 0));
}

TEST(ErrorHandlingScalar, CheckNonNegativeVectorization) {
  Eigen::MatrixXd m = Eigen::MatrixXd::Constant(3, 2, 0);
  EXPECT_NO_THROW(
      check_nonnegative(function, "m", std::vector<Eigen::MatrixXd>{m, m, m}));
  Eigen::MatrixXd m2 = m;
  m2(1, 1) = -1;
  EXPECT_THROW(
      check_nonnegative(function, "m", std::vector<Eigen::MatrixXd>{m, m2, m}),
      std::domain_error);
}
