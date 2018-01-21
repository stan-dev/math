#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

using stan::math::check_less_or_equal;

TEST(ErrorHandlingScalar, CheckLessOrEqual_Matrix) {
  const char* function = "check_less_or_equal";
  double x;
  double high;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<double, Eigen::Dynamic, 1> high_vec;
  x_vec.resize(3);
  high_vec.resize(3);

  // x_vec, high
  x_vec << -5, 0, 5;
  high = 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, 5;
  high = 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = 5;
  EXPECT_THROW(check_less_or_equal(function, "x", x_vec, high),
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high));

  // x_vec, high_vec
  x_vec << -5, 0, 5;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << std::numeric_limits<double>::infinity(), 10, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, 5;
  high_vec << 10, 10, 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, 10;
  EXPECT_THROW(check_less_or_equal(function, "x", x_vec, high_vec),
               std::domain_error);

  x_vec << -5, 0, std::numeric_limits<double>::infinity();
  high_vec << 10, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x_vec, high_vec));

  // x, high_vec
  x = -100;
  high_vec << 0, 5, 10;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = 10;
  high_vec << 100, 200, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = 5;
  high_vec << 100, 200, 5;
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));

  x = std::numeric_limits<double>::infinity();
  high_vec << 10, 20, 30;
  EXPECT_THROW(check_less_or_equal(function, "x", x, high_vec),
               std::domain_error);

  x = std::numeric_limits<double>::infinity();
  high_vec << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less_or_equal(function, "x", x, high_vec));
}

TEST(ErrorHandlingScalar, CheckLessOrEqual_Matrix_one_indexed_message) {
  const char* function = "check_less";
  double x;
  double high;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<double, Eigen::Dynamic, 1> high_vec;
  x_vec.resize(3);
  high_vec.resize(3);
  std::string message;

  // x_vec, high
  x_vec << -5, 0, 5;
  high = 4;

  try {
    check_less_or_equal(function, "x", x_vec, high);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]")) << message;

  // x_vec, high_vec
  x_vec << -5, 6, 0;
  high_vec << 10, 5, 10;

  try {
    check_less_or_equal(function, "x", x_vec, high_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[2]")) << message;

  // x, high_vec
  x = 30;
  high_vec << 10, 20, 30;

  try {
    check_less_or_equal(function, "x", x, high_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_EQ(std::string::npos, message.find("["))
      << "no index provided" << std::endl
      << message;
}

TEST(ErrorHandlingScalar, CheckLessOrEqual_nan) {
  const char* function = "check_less_or_equal";
  double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec(3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> low_vec(3);

  // x_vec, low_vec
  x_vec << -1, 0, 1;
  low_vec << -2, -1, 0;
  EXPECT_THROW(check_less_or_equal(function, "x", x_vec, nan),
               std::domain_error);

  for (int i = 0; i < x_vec.size(); i++) {
    x_vec << -1, 0, 1;
    x_vec(i) = nan;
    EXPECT_THROW(check_less_or_equal(function, "x", x_vec, low_vec),
                 std::domain_error);

    x_vec << -1, 0, 1;
    for (int i = 0; i < low_vec.size(); i++) {
      low_vec << -1, 0, 1;
      low_vec(i) = nan;
      EXPECT_THROW(check_less_or_equal(function, "x", x_vec, low_vec),
                   std::domain_error);
    }

    for (int i = 0; i < x_vec.size(); i++) {
      x_vec << -1, 0, 1;
      low_vec << -2, -1, 0;
      x_vec(i) = nan;
      for (int j = 0; j < low_vec.size(); j++) {
        low_vec(i) = nan;
        EXPECT_THROW(check_less_or_equal(function, "x", x_vec, low_vec),
                     std::domain_error);
      }
    }
  }
}
