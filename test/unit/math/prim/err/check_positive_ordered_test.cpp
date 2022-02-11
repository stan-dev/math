#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, checkPositiveOrdered) {
  using stan::math::check_positive_ordered;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(3);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << 0, 1, 2;
  }
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));

  for (auto& y_i : y) {
    y_i << 0, 10, std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));

  for (auto& y_i : y) {
    y_i << 0, 0, 0;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << 0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << -1, 0, 0;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << 0, 3, 2;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkPositiveOrdered_one_indexed_message) {
  using stan::math::check_positive_ordered;
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(3);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << -1, 0, 0;
  }
  try {
    check_positive_ordered("check_positive_ordered", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("element at 1")) << message;

  for (auto& y_i : y) {
    y_i << 0, 5, 1;
  }
  try {
    check_positive_ordered("check_positive_ordered", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("element at 3")) << message;
}

TEST(ErrorHandlingMatrix, checkPositiveOrdered_nan) {
  using stan::math::check_positive_ordered;
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(3);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};

  for (auto& y_i : y) {
    y_i << 0, 1, 2;
  }
  for (int i = 0; i < y.size(); i++) {
    y[i][i] = nan;
    EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
                 std::domain_error);
    y[i][i] = i;
  }
  for (int i = 0; i < y.size(); i++) {
    y[i] << 0, 10, std::numeric_limits<double>::infinity();
    y[i][i] = nan;
    EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
                 std::domain_error);
    y[i][i] = i;
  }
}
