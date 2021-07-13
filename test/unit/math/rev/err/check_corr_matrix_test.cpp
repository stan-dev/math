#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, CheckCorrVarMatrix) {
  using stan::math::check_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  using var_mat = stan::math::var_value<Eigen::MatrixXd>;
  y.resize(2, 2);

  y << 1, 0, 0, 1;
  EXPECT_NO_THROW(check_corr_matrix("test", "y", var_mat(y)));

  y << 10, 0, 0, 10;
  EXPECT_THROW(check_corr_matrix("test", "y", var_mat(y)), std::domain_error);
}

TEST(ErrorHandlingMatrix, CheckCorrVarMatrix_one_indexed_message) {
  using stan::math::check_corr_matrix;
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  using var_mat = stan::math::var_value<Eigen::MatrixXd>;
  y.resize(2, 2);

  y << 10, 0, 0, 1;
  try {
    check_corr_matrix("test", "y", var_mat(y));
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("(1,1)")) << message;

  EXPECT_EQ(std::string::npos, message.find("(0, 0)")) << message;
}

TEST(ErrorHandlingMatrix, CheckCorrVarMatrix_nan) {
  using stan::math::check_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y;
  using var_mat = stan::math::var_value<Eigen::MatrixXd>;
  y.resize(2, 2);
  double nan = std::numeric_limits<double>::quiet_NaN();

  for (int i = 0; i < y.size(); i++) {
    y << 1, 0, 0, 1;
    y(i) = nan;
    EXPECT_THROW(check_corr_matrix("test", "y", var_mat(y)), std::domain_error);

    y << 10, 0, 0, 10;
    y(i) = nan;
    EXPECT_THROW(check_corr_matrix("test", "y", var_mat(y)), std::domain_error);
  }
}
