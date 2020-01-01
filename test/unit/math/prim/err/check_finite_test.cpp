#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

using stan::math::check_finite;

// ---------- check_finite: vector tests ----------
TEST(ErrorHandlingArr, CheckFinite_Vector) {
  const char* function = "check_finite";
  std::vector<double> x = {-1, 0, 1};
  ASSERT_NO_THROW(check_finite(function, "x", x))
      << "check_finite should be true with finite x";

  x = {-1, 0, std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on Inf";

  x = {-1, 0, -std::numeric_limits<double>::infinity()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on -Inf";

  x = {-1, 0, std::numeric_limits<double>::quiet_NaN()};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingArr, CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {nan, 0, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, nan, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, 0, nan};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);
}

// ---------- check_finite: matrix tests ----------
TEST(ErrorHandlingMat, CheckFinite_Matrix) {
  const char* function = "check_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;
  ASSERT_NO_THROW(check_finite(function, "x", x))
      << "check_finite should be true with finite x";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on Inf";

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on -Inf";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingMat, CheckFinite_Matrix_one_indexed_message) {
  const char* function = "check_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  std::string message;

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  try {
    check_finite(function, "x", x);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[3]")) << message;
}

TEST(ErrorHandlingMat, CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_mat(3);
  x_mat << nan, 0, 1;
  EXPECT_THROW(check_finite(function, "x_mat", x_mat), std::domain_error);

  x_mat << 1, nan, 1;
  EXPECT_THROW(check_finite(function, "x_mat", x_mat), std::domain_error);

  x_mat << 1, 0, nan;
  EXPECT_THROW(check_finite(function, "x_mat", x_mat), std::domain_error);
}

TEST(ErrorHandlingScalar, CheckFinite) {
  const char* function = "check_finite";
  double x = 0;

  EXPECT_NO_THROW(check_finite(function, "x", x))
      << "check_finite should be true with finite x: " << x;
  x = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on Inf: " << x;
  x = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on -Inf: " << x;

  x = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error)
      << "check_finite should throw exception on NaN: " << x;
}

TEST(ErrorHandlingScalar, CheckFinite_nan) {
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_finite(function, "x", nan), std::domain_error);
}
