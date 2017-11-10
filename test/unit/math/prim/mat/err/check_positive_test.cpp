#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingScalar,CheckPositive) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  Eigen::Matrix<double,Eigen::Dynamic,1> x_mat(3);
  x_mat << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    EXPECT_NO_THROW(check_positive(function, "x", x_mat));
  }

  x_mat(0) = 0;
  EXPECT_THROW(check_positive(function, "x", x_mat),
               std::domain_error);
}

TEST(ErrorHandlingScalar,CheckPositive_nan) {
  using stan::math::check_positive;
  const char* function = "check_positive";

  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double,Eigen::Dynamic,1> x_mat(3);
  x_mat   << 1, 2, 3;
  for (int i = 0; i < x_mat.size(); i++) {
    x_mat(i) = nan;
    EXPECT_THROW(check_positive(function, "x", x_mat),
                 std::domain_error);
    x_mat(i) = i;
  }
}
