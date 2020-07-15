#include <stan/math/prim/err/check_sorted.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

using stan::math::check_sorted;

TEST(ErrorHandling, checkSorted) {
  std::vector<double> y = {0, 1, 2};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {0, 1, std::numeric_limits<double>::infinity()};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {-10, 1, std::numeric_limits<double>::infinity()};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {-std::numeric_limits<double>::infinity(), 1,
       std::numeric_limits<double>::infinity()};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {0, 0, 0};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {0, std::numeric_limits<double>::infinity(),
       std::numeric_limits<double>::infinity()};
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y = {-1, 3, 2};
  EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
}

TEST(ErrorHandling, checkSorted_one_indexed_message) {
  std::string message;
  std::vector<double> y = {0, 5, 1};

  try {
    check_sorted("check_sorted", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("element at 3")) << message;
}

TEST(ErrorHandling, checkSorted_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> y = {0, 1, 2};

  for (size_t i = 0; i < y.size(); i++) {
    y[i] = nan;
    EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
    y[i] = i;
  }
  for (size_t i = 0; i < y.size(); i++) {
    y = {0.0, 10.0, std::numeric_limits<double>::infinity()};
    y[i] = nan;
    EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
  }
}

TEST(ErrorHandlingMatrix, checkSorted) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);

  y << 0, 1, 2;
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << 0, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << -10, 10, std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << -std::numeric_limits<double>::infinity(), 10,
      std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << 0, 0, 0;
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << 0, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_sorted("check_sorted", "y", y));

  y << -1, 3, 2;
  EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
}

TEST(ErrorHandlingMatrix, checkSorted_one_indexed_message) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  y.resize(3);

  y << 0, 5, 1;
  try {
    check_sorted("check_sorted", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("element at 3")) << message;
}

TEST(ErrorHandlingMatrix, checkSorted_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  double nan = std::numeric_limits<double>::quiet_NaN();
  y.resize(3);

  y << 0, 1, 2;
  for (int i = 0; i < y.size(); i++) {
    y[i] = nan;
    EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
    y[i] = i;
  }
  for (int i = 0; i < y.size(); i++) {
    y << 0, 10, std::numeric_limits<double>::infinity();
    y[i] = nan;
    EXPECT_THROW(check_sorted("check_sorted", "y", y), std::domain_error);
  }
}
