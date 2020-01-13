#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ErrorHandlingArr, checkNonzeroSizeMatrix) {
  using stan::math::check_nonzero_size;

  std::vector<double> a(4, 3);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "a", a));

  a.resize(2);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "a", a));

  a.resize(0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "a", a),
                   std::invalid_argument, "has size 0");
}

TEST(ErrorHandlingArr, checkNonzeroSizeMatrix_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> a(4, nan);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "a", a));

  a.resize(2);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "a", a));

  a.resize(0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "a", a),
                   std::invalid_argument, "has size 0");
}

TEST(ErrorHandlingMatrix, checkNonzeroSizeMatrix) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  using stan::math::check_nonzero_size;

  y.resize(3, 3);
  EXPECT_NO_THROW(check_nonzero_size("checkNonzeroSize", "y", y));
  y.resize(2, 3);
  EXPECT_NO_THROW(check_nonzero_size("checkNonzeroSize", "y", y));

  y.resize(0, 0);
  EXPECT_THROW_MSG(check_nonzero_size("checkNonzeroSize", "y", y),
                   std::invalid_argument, "has size 0");
}

TEST(ErrorHandlingMatrix, checkNonzeroSizeMatrix_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();

  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "y", y));
  y.resize(2, 3);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize", "y", y));

  y.resize(0, 0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "y", y),
                   std::invalid_argument, "has size 0");
}
