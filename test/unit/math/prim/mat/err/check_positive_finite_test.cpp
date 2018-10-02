#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

using stan::math::check_positive_finite;

TEST(ErrorHandlingScalar, CheckPositiveFinite_Matrix) {
  const char* function = "check_positive_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << 3, 2, 1;
  ASSERT_NO_THROW(check_positive_finite(function, "x", x))
      << "check_positive_finite should be true with finite x";

  x.resize(3);
  x << 2, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on Inf";

  x.resize(3);
  x << 0, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on x=0";

  x.resize(3);
  x << -1, 1, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on x=-1";

  x.resize(3);
  x << 2, 1, -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on -Inf";

  x.resize(3);
  x << 1, 2, std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_positive_finite(function, "x", x), std::domain_error)
      << "check_positive_finite should throw exception on NaN";
}

TEST(ErrorHandlingScalar, CheckPositiveFinite_Matrix_one_indexed_message) {
  const char* function = "check_positive_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  std::string message;

  x.resize(3);
  x << 1, 2, std::numeric_limits<double>::infinity();
  try {
    check_positive_finite(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]")) << message;
}
TEST(ErrorHandlingScalar, CheckPositiveFinite_Matrix_one_indexed_message_2) {
  const char* function = "check_positive_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  std::string message;

  x.resize(3);
  x << -1, 2, std::numeric_limits<double>::infinity();
  try {
    check_positive_finite(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[1]")) << message;
}

TEST(ErrorHandlingScalar, CheckPositiveFinite_Matrix_one_indexed_message_3) {
  const char* function = "check_positive_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  std::string message;

  x.resize(3);
  x << 1, 0, std::numeric_limits<double>::infinity();
  try {
    check_positive_finite(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[2]")) << message;
}

TEST(ErrorHandlingScalar, CheckPositiveFinite_nan) {
  const char* function = "check_positive_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    x_mat(i) = nan;
    EXPECT_THROW(check_positive_finite(function, "x", x_mat),
                 std::domain_error);
    x_mat(i) = i;
  }
}
