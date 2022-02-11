#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMat, CheckGreater_Matrix) {
  using stan::math::check_greater;
  const char* function = "check_greater";
  double x;
  double low;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_mat;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> low_mat;
  x_mat.resize(3, 3);
  low_mat.resize(3, 3);
  std::vector<double> x_scalar_vec{x, x, x};
  std::vector<double> low_scalar_vec{low, low, low};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x_vec{
      x_mat, x_mat, x_mat};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> low_vec{
      low_mat, low_mat, low_mat};
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i] << -1, 0, 1, -1, 0, 1, -1, 0, 1;
    low_vec[i] << -2, -1, 0, -2, -1, 0, -2, -1, 0;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>";

  for (int i = 0; i < x_vec.size(); ++i) {
    low_vec[i] << -1.1, -0.1, 0.9, -1.1, -0.1, 0.9, -1.1, -0.1, 0.9;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, std::numeric_limits<double>::infinity();
    low_vec[i].col(i) << -2, -1, 0;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>, y has infinity";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_vec[i].col(i) << -2, 0, 0;
  }
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should fail for index 1";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_vec[i].col(i) << -2, -1, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should fail with infinity";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, std::numeric_limits<double>::infinity();
    low_vec[i].col(i) << -2, -1, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_greater(function, "x", x_vec, low_vec), std::domain_error)
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should fail with infinity";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_vec[i].col(i) << -2, -1, -std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_vec))
      << "check_greater: matrix<3, 1>, matrix<3, 1>, should pass with "
         "-infinity";

  // x_vec, low
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_scalar_vec[i] = -2;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_scalar_vec))
      << "check_greater: matrix<3, 1>, double";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_scalar_vec[i] = -std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_vec, low_scalar_vec))
      << "check_greater: matrix<3, 1>, double";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_scalar_vec[i] = 0;
  }
  EXPECT_THROW(check_greater(function, "x", x_vec, low_scalar_vec),
               std::domain_error)
      << "check_greater: matrix<3, 1>, double, should fail for index 1/2";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -1, 0, 1;
    low_scalar_vec[i] = std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_greater(function, "x", x, low_scalar_vec),
               std::domain_error)
      << "check_greater: matrix<3, 1>, double, should fail with infinity";

  // x, low_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 2;
    low_vec[i].col(i) << -1, 0, 1;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_scalar_vec, low_vec))
      << "check_greater: double, matrix<3, 1>";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 10;
    low_vec[i].col(i) << -1, 0, -std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_scalar_vec, low_vec))
      << "check_greater: double, matrix<3, 1>, low has -inf";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 10;
    low_vec[i].col(i) << -1, 0, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_greater(function, "x", x_scalar_vec, low_vec),
               std::domain_error)
      << "check_greater: double, matrix<3, 1>, low has inf";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = std::numeric_limits<double>::infinity();
    low_vec[i].col(i) << -1, 0, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_greater(function, "x", x_scalar_vec, low_vec),
               std::domain_error)
      << "check_greater: double, matrix<3, 1>, x is inf, low has inf";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = std::numeric_limits<double>::infinity();
    low_vec[i].col(i) << -1, 0, 1;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_scalar_vec, low_vec))
      << "check_greater: double, matrix<3, 1>, x is inf";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 1.1;
    low_vec[i].col(i) << -1, 0, 1;
  }
  EXPECT_NO_THROW(check_greater(function, "x", x_scalar_vec, low_vec))
      << "check_greater: double, matrix<3, 1>";

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 0.9;
    low_vec[i].col(i) << -1, 0, 1;
  }
  EXPECT_THROW(check_greater(function, "x", x_scalar_vec, low_vec),
               std::domain_error)
      << "check_greater: double, matrix<3, 1>";
}

TEST(ErrorHandlingMat, CheckGreater_Matrix_one_indexed_message) {
  using stan::math::check_greater;
  static const char* function = "check_greater";
  double x;
  double low;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_mat;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> low_mat;
  x_mat.resize(3, 3);
  low_mat.resize(3, 3);
  std::vector<double> x_scalar_vec{x, x, x};
  std::vector<double> low_scalar_vec{low, low, low};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x_vec{
      x_mat, x_mat, x_mat};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> low_vec{
      low_mat, low_mat, low_mat};

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].resize(3, 3);
    low_vec[i].resize(3, 3);
  }
  std::string message;

  // x_vec, low
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i] << 10, 10, -1, 10, 10, -1, 10, 10, -1;
    low_scalar_vec[i] = 0;
  }

  try {
    check_greater(function, "x", x_vec, low_scalar_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[1][1, 3]")) << message;

  // x_vec, low_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i] << 10, -1, 10, 10, -1, 10, 10, -1, 10;
    low_vec[i] << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  }

  try {
    check_greater(function, "x", x_vec, low_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[1][1, 2]")) << message;

  // x, low_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 0;
    low_vec[i] << -1, 10, 10, -1, 10, 10, -1, 10, 10;
  }

  try {
    check_greater(function, "x", x_scalar_vec, low_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[1]"))
      << "no index provided" << std::endl
      << message;
}

TEST(ErrorHandlingMat, CheckGreater_nan) {
  using stan::math::check_greater;
  const char* function = "check_greater";
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_greater(function, "x", nan, lb), std::domain_error);
  EXPECT_THROW(check_greater(function, "x", x, nan), std::domain_error);
  EXPECT_THROW(check_greater(function, "x", nan, nan), std::domain_error);
}

TEST(ErrorHandlingScalar, CheckGreater) {
  using stan::math::check_greater;
  const char* function = "check_greater";
  double x = 10.0;
  double lb = 0.0;

  EXPECT_NO_THROW(check_greater(function, "x", x, lb))
      << "check_greater should be true with x > lb";

  x = -1.0;
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x < lb";

  x = lb;
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == lb";

  x = std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_greater(function, "x", x, lb))
      << "check_greater should be true with x == Inf and lb = 0.0";

  x = 10.0;
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == 10.0 and lb == Inf";

  x = std::numeric_limits<double>::infinity();
  lb = std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_greater(function, "x", x, lb), std::domain_error)
      << "check_greater should throw an exception with x == Inf and lb == Inf";
}

TEST(ErrorHandlingScalar, CheckGreater_nan) {
  using stan::math::check_greater;
  const char* function = "check_greater";
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_greater(function, "x", nan, lb), std::domain_error);
  EXPECT_THROW(check_greater(function, "x", x, nan), std::domain_error);
  EXPECT_THROW(check_greater(function, "x", nan, nan), std::domain_error);
}
