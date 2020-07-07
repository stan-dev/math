#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMat, CheckGreaterOrEqualMatrix) {
  using stan::math::check_greater_or_equal;
  const char* function = "check_greater_or_equal";
  double x;
  double low;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<double, Eigen::Dynamic, 1> low_vec;
  x_vec.resize(3);
  low_vec.resize(3);

  // x_vec, low_vec
  x_vec << -1, 0, 1;
  low_vec << -2, -1, 0;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>";

  x_vec << -1, 0, 1;
  low_vec << -1.1, -0.1, 0.9;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>";

  x_vec << -1, 0, std::numeric_limits<double>::infinity();
  low_vec << -2, -1, 0;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>, y has infinity";

  x_vec << -1, 0, 1;
  low_vec << -2, 0, 0;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>, "
      << "should pass for index 1";

  x_vec << -1, 0, 1;
  low_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low_vec),
               std::domain_error)
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>, "
      << "should fail with infinity";

  x_vec << -1, 0, std::numeric_limits<double>::infinity();
  low_vec << -2, -1, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>, both bound and "
      << "value infinity";

  x_vec << -1, 0, 1;
  low_vec << -2, -1, -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low_vec))
      << "check_greater_or_equal: matrix<3, 1>, matrix<3, 1>, should pass with "
      << "-infinity";

  // x_vec, low
  x_vec << -1, 0, 1;
  low = -2;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low))
      << "check_greater_or_equal: matrix<3, 1>, double";

  x_vec << -1, 0, 1;
  low = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x_vec, low))
      << "check_greater_or_equal: matrix<3, 1>, double";

  x_vec << -1, 0, 1;
  low = 0;
  EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low),
               std::domain_error)
      << "check_greater_or_equal: matrix<3, 1>, double, should fail for "
         "index 1/2";

  x_vec << -1, 0, 1;
  low = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low),
               std::domain_error)
      << "check_greater_or_equal: matrix<3, 1>, double, should fail with "
         "infinity";

  // x, low_vec
  x = 2;
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, low_vec))
      << "check_greater_or_equal: double, matrix<3, 1>";

  x = 10;
  low_vec << -1, 0, -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, low_vec))
      << "check_greater_or_equal: double, matrix<3, 1>, low has -inf";

  x = 10;
  low_vec << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x, low_vec),
               std::domain_error)
      << "check_greater_or_equal: double, matrix<3, 1>, low has inf";

  x = std::numeric_limits<double>::infinity();
  low_vec << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, low_vec))
      << "check_greater_or_equal: double, matrix<3, 1>, x is inf, low has inf";

  x = std::numeric_limits<double>::infinity();
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, low_vec))
      << "check_greater_or_equal: double, matrix<3, 1>, x is inf";

  x = 1.1;
  low_vec << -1, 0, 1;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, low_vec))
      << "check_greater_or_equal: double, matrix<3, 1>";

  x = 0.9;
  low_vec << -1, 0, 1;
  EXPECT_THROW(check_greater_or_equal(function, "x", x, low_vec),
               std::domain_error)
      << "check_greater_or_equal: double, matrix<3, 1>";
}

TEST(ErrorHandlingMat, CheckGreaterOrEqual_Matrix_one_indexed_message) {
  using stan::math::check_greater_or_equal;
  const char* function = "check_greater_or_equal";
  double x;
  double low;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec;
  Eigen::Matrix<double, Eigen::Dynamic, 1> low_vec;
  x_vec.resize(3);
  low_vec.resize(3);
  std::string message;

  // x_vec, low
  x_vec << 10, 10, -1;
  low = 0;

  try {
    check_greater_or_equal(function, "x", x_vec, low);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]")) << message;

  // x_vec, low_vec
  x_vec << 10, -1, 10;
  low_vec << 0, 0, 0;

  try {
    check_greater_or_equal(function, "x", x_vec, low_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[2]")) << message;

  // x, low_vec
  x = 0;
  low_vec << -1, 10, 10;

  try {
    check_greater_or_equal(function, "x", x, low_vec);
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

TEST(ErrorHandlingMat, CheckGreaterOrEqual_nan) {
  using stan::math::check_greater_or_equal;
  const char* function = "check_greater_or_equal";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec(3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> low_vec(3);

  // x_vec, low_vec
  x_vec << -1, 0, 1;
  low_vec << -2, -1, 0;
  EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, nan),
               std::domain_error);

  for (int i = 0; i < x_vec.size(); i++) {
    x_vec << -1, 0, 1;
    x_vec(i) = nan;
    EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low_vec),
                 std::domain_error);
  }

  x_vec << -1, 0, 1;
  for (int i = 0; i < low_vec.size(); i++) {
    low_vec << -1, 0, 1;
    low_vec(i) = nan;
    EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low_vec),
                 std::domain_error);
  }

  for (int i = 0; i < x_vec.size(); i++) {
    x_vec << -1, 0, 1;
    low_vec << -2, -1, 0;
    x_vec(i) = nan;
    for (int j = 0; j < low_vec.size(); j++) {
      low_vec(i) = nan;
      EXPECT_THROW(check_greater_or_equal(function, "x", x_vec, low_vec),
                   std::domain_error);
    }
  }
}

TEST(ErrorHandlingScalar, CheckGreaterOrEqual) {
  using stan::math::check_greater_or_equal;
  const char* function = "check_greater_or_equal";
  double x = 10.0;
  double lb = 0.0;

  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
      << "check_greater_or_equal should be true with x > lb";

  x = -1.0;
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb), std::domain_error)
      << "check_greater_or_equal should throw an exception with x < lb";

  x = lb;
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
      << "check_greater_or_equal should not throw an exception with x == lb";

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
      << "check_greater should be true with x == Inf and lb = 0.0";

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater_or_equal(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == 10.0 and lb == Inf";

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater_or_equal(function, "x", x, lb))
      << "check_greater should not throw an exception with "
      << "x == Inf and lb == Inf";
}

TEST(ErrorHandlingScalar, CheckGreaterOrEqual_nan) {
  using stan::math::check_greater_or_equal;
  const char* function = "check_greater_or_equal";
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_greater_or_equal(function, "x", nan, lb),
               std::domain_error);
  EXPECT_THROW(check_greater_or_equal(function, "x", x, nan),
               std::domain_error);
  EXPECT_THROW(check_greater_or_equal(function, "x", nan, nan),
               std::domain_error);
}
