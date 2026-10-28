#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, checkPositiveOrdered) {
  using stan::math::check_positive_ordered;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(3);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << 0.5, 1, 2;
  }
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));

  for (auto& y_i : y) {
    y_i << 1, 10, std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));

  for (auto& y_i : y) {
    y_i << std::numeric_limits<double>::denorm_min(), 1, 2;
  }
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));

  for (auto& y_i : y) {
    y_i << 0, 1, 2;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << -0.0, 1, 2;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << 0, 10, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << 1, 1, 1;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i << 1, std::numeric_limits<double>::infinity(),
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
    y_i << 1, 3, 2;
  }
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
               std::domain_error);

  Eigen::Matrix<double, Eigen::Dynamic, 1> y_one(1);
  y_one << 1;
  EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y_one));
  y_one << 0;
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y_one),
               std::domain_error);
  y_one << std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y_one),
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
    y_i << 0, 1, 2;
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
  EXPECT_NE(std::string::npos, message.find("should be positive")) << message;

  for (auto& y_i : y) {
    y_i << 1, 5, 2;
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
  Eigen::Matrix<double, Eigen::Dynamic, 1> base1(3);
  base1 << 1, 2, 3;
  Eigen::Matrix<double, Eigen::Dynamic, 1> base2(3);
  base2 << 1, 10, std::numeric_limits<double>::infinity();

  for (const auto& base : {base1, base2}) {
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{base, base, base};
    EXPECT_NO_THROW(check_positive_ordered("check_positive_ordered", "y", y));
    for (int i = 0; i < y.size(); i++) {
      y[i][i] = nan;
      EXPECT_THROW(check_positive_ordered("check_positive_ordered", "y", y),
                   std::domain_error);
      y[i] = base;
    }
  }
}
