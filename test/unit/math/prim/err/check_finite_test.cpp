#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(ErrorHandlingArr, CheckFinite_Vector) {
  using stan::math::check_finite;
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

TEST(ErrorHandlingArr, CheckFinite_std_vector_std_vector) {
  using stan::math::check_finite;
  const char* function = "check_finite";
  std::vector<double> x = {-1, 0, 1};
  std::vector<std::vector<double>> xx = {x};
  ASSERT_NO_THROW(check_finite(function, "x", xx))
      << "check_finite should be true with finite x";

  x = {-1, 0, std::numeric_limits<double>::infinity()};
  xx = {x};
  EXPECT_THROW(check_finite(function, "x", xx), std::domain_error)
      << "check_finite should throw exception on Inf";

  x = {-1, 0, -std::numeric_limits<double>::infinity()};
  xx = {x};
  EXPECT_THROW(check_finite(function, "x", xx), std::domain_error)
      << "check_finite should throw exception on -Inf";

  x = {-1, 0, std::numeric_limits<double>::quiet_NaN()};
  xx = {x};
  EXPECT_THROW(check_finite(function, "x", xx), std::domain_error)
      << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingArr, CheckFinite_nan) {
  using stan::math::check_finite;
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x = {nan, 0, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, nan, 1};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);

  x = {1, 0, nan};
  EXPECT_THROW(check_finite(function, "x", x), std::domain_error);
}

TEST(ErrorHandlingMat, CheckFinite_Matrix) {
  using stan::math::check_finite;
  const char* function = "check_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;
  ASSERT_NO_THROW(check_finite(function, "x", x))
      << "check_finite should be true with finite x";

  ASSERT_NO_THROW(check_finite(function, "x", x.array()))
      << "check_finite should be true with finite x";

  ASSERT_NO_THROW(check_finite(function, "x", x.transpose()))
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

TEST(ErrorHandlingMat, CheckFinite_std_vector_Matrix) {
  using stan::math::check_finite;
  const char* function = "check_finite";
  Eigen::Matrix<double, Eigen::Dynamic, 1> x;

  x.resize(3);
  x << -1, 0, 1;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> xv = {x};
  std::vector<Eigen::Array<double, Eigen::Dynamic, 1>> xva = {x.array()};
  std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic>> xvt = {x.transpose()};

  ASSERT_NO_THROW(check_finite(function, "x", xv))
      << "check_finite should be true with finite x";

  ASSERT_NO_THROW(check_finite(function, "x", xva))
      << "check_finite should be true with finite x";

  ASSERT_NO_THROW(check_finite(function, "x", xvt))
      << "check_finite should be true with finite x";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::infinity();
  xv = {x};
  EXPECT_THROW(check_finite(function, "x", xv), std::domain_error)
      << "check_finite should throw exception on Inf";

  x.resize(3);
  x << -1, 0, -std::numeric_limits<double>::infinity();
  xv = {x};
  EXPECT_THROW(check_finite(function, "x", xv), std::domain_error)
      << "check_finite should throw exception on -Inf";

  x.resize(3);
  x << -1, 0, std::numeric_limits<double>::quiet_NaN();
  xv = {x};
  EXPECT_THROW(check_finite(function, "x", xv), std::domain_error)
      << "check_finite should throw exception on NaN";
}

TEST(ErrorHandlingMat, CheckFinite_Matrix_one_indexed_message) {
  using stan::math::check_finite;
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
  using stan::math::check_finite;
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
  using stan::math::check_finite;
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
  using stan::math::check_finite;
  const char* function = "check_finite";
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_finite(function, "x", nan), std::domain_error);
}
