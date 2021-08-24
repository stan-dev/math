#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, CheckCorrMatrix) {
  using stan::math::check_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};

  for (auto& y : y_vec) {
    y.resize(2, 2);
    y << 1, 0, 0, 1;
  }
  EXPECT_NO_THROW(check_corr_matrix("test", "y", y_vec));

  for (auto& y : y_vec) {
    y << 10, 0, 0, 10;
  }
  EXPECT_THROW(check_corr_matrix("test", "y", y_vec), std::domain_error);
}

TEST(ErrorHandlingMatrix, CheckCorrMatrix_one_indexed_message) {
  using stan::math::check_corr_matrix;
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};
  for (auto& y : y_vec) {
    y.resize(2, 2);
    y << 10, 0, 0, 1;
  }
  try {
    check_corr_matrix("test", "y", y_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("(1,1)")) << message;

  EXPECT_EQ(std::string::npos, message.find("(0, 0)")) << message;
}

TEST(ErrorHandlingMatrix, CheckCorrMatrix_nan) {
  using stan::math::check_corr_matrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_mat;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> y_vec{
      y_mat, y_mat, y_mat};
  for (auto& y : y_vec) {
    y.resize(2, 2);
  }
  double nan = std::numeric_limits<double>::quiet_NaN();

  for (int i = 0; i < y_vec[0].size(); i++) {
    for (auto& y : y_vec) {
      y << 1, 0, 0, 1;
      y(i) = nan;
    }
    EXPECT_THROW(check_corr_matrix("test", "y", y_vec), std::domain_error);

    for (auto& y : y_vec) {
      y << 10, 0, 0, 10;
      y(i) = nan;
    }
    EXPECT_THROW(check_corr_matrix("test", "y", y_vec), std::domain_error);
  }
}
