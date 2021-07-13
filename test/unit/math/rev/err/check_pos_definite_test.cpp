#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(AgradRevErrorHandlingMatrix, checkPosDefiniteMatrix_nan) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;

  Matrix<var, Dynamic, Dynamic> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  using stan::math::check_pos_definite;

  y.resize(1, 1);
  y << nan;
  EXPECT_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y),
               std::domain_error);

  y.resize(3, 3);
  y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  EXPECT_NO_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y));

  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++) {
      y << 2, -1, 0, -1, 2, -1, 0, -1, 2;
      y(i, j) = nan;
      if (i >= j)
        EXPECT_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y),
                     std::domain_error);
    }

  y << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  EXPECT_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y),
               std::domain_error);
}

TEST(AgradRevErrorHandlingMatrix, checkPosDefiniteMatrix_nan_varmat) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var_value;

  Matrix<double, Dynamic, Dynamic> y_val(3, 3);
  y_val << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  var_value<Matrix<double, Dynamic, Dynamic>> y(y_val);
  using stan::math::check_pos_definite;

  EXPECT_NO_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y));

  double nan = std::numeric_limits<double>::quiet_NaN();
  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++) {
      y.vi_->val_(i, j) = nan;
      if (i >= j)
        EXPECT_THROW(check_pos_definite("checkPosDefiniteMatrix", "y", y),
                     std::domain_error);
      y.vi_->val_ = y_val;
    }

  Matrix<double, Dynamic, Dynamic> y_val_zero(3, 3);
  y_val_zero << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  var_value<Matrix<double, Dynamic, Dynamic>> y_zero(y_val_zero);
  EXPECT_THROW(check_pos_definite("checkPosDefiniteMatrix", "y_zero", y_zero),
               std::domain_error);
}
