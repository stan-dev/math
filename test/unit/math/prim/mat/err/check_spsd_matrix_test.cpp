#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(ErrorHandlingMatrix, checkSpsdMatrixPosDef) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_NO_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix",
                                                "y", y));

  y << 1, 2, 3, 2, 1, 2, 3, 2, 1;
  EXPECT_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix", "y", y),
               std::domain_error);

  y.setZero();
  EXPECT_NO_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix", "y", y));
}

TEST(ErrorHandlingMatrix, checkSpsdMatrixZero) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 3);
  EXPECT_NO_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix", "y", y));
}

TEST(ErrorHandlingMatrix, checkSpsdNotSquare) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3, 2);
  EXPECT_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix", "y", y),
               std::invalid_argument);
}

TEST(ErrorHandlingMatrix, checkSpsdMatrixPosDef_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_NO_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix",
                                                "y", y));

  y.setZero();
  EXPECT_NO_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix", "y", y));

  for (int i = 0; i < y.size(); i++) {
    y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
    y(i) = nan;
    EXPECT_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix",
                                               "y", y),
                 std::domain_error);

    y.setZero();
    y(i) = nan;
    EXPECT_THROW(stan::math::check_spsd_matrix("checkSpsdMatrix",
                                               "y", y),
                 std::domain_error);
  }
}
