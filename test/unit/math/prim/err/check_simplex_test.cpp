#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <string>

TEST(ErrorHandlingMatrix, checkSimplex) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i << 0.5, 0.5;
  }

  EXPECT_NO_THROW(stan::math::check_simplex("checkSimplex", "y", y));

  for (auto& y_i : y) {
    y_i[1] = 0.55;
  }
  EXPECT_THROW(stan::math::check_simplex("checkSimplex", "y", y),
               std::domain_error);
}

TEST(ErrorHandlingMatrix, checkSimplex_message_negative_value) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(100);
  y_vec.setZero();
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i[0] = -0.1;
    y_i[1] = 1.1;
  }

  try {
    stan::math::check_simplex("checkSimplex", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos != message.find(" y[1] is not a valid simplex"))
      << message;

  EXPECT_TRUE(std::string::npos != message.find("y[1][1] = -0.1")) << message;

  for (auto& y_i : y) {
    y_i.setZero();
    y_i[0] = 0.1;
    y_i[1] = -0.1;
    y_i[2] = 1.0;
  }
  try {
    stan::math::check_simplex("checkSimplex", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos != message.find(" y[1] is not a valid simplex"))
      << message;

  EXPECT_TRUE(std::string::npos != message.find(" y[1][2] = -0.1")) << message;
}

TEST(ErrorHandlingMatrix, checkSimplex_message_sum) {
  std::string message;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(100);
  y_vec.setZero();
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  for (auto& y_i : y) {
    y_i[13] = 0.9;
  }

  try {
    stan::math::check_simplex("checkSimplex", "y", y);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_TRUE(std::string::npos != message.find(" y[1] is not a valid simplex"))
      << message;

  EXPECT_TRUE(std::string::npos != message.find("sum(y[1]) = 0.9")) << message;
}

TEST(ErrorHandlingMatrix, checkSimplex_message_length) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(0);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};

  using stan::math::check_simplex;

  EXPECT_THROW_MSG(check_simplex("checkSimplex", "y", y), std::invalid_argument,
                   "y[1] has size 0, but must have a non-zero size");
}

TEST(ErrorHandlingMatrix, checkSimplex_nan) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> y_vec(2);
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> y{y_vec, y_vec, y_vec};
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  for (auto& y_i : y) {
    y_i << nan, 0.5;
  }

  EXPECT_THROW(stan::math::check_simplex("checkSimplex", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[1] = 0.55;
  }
  EXPECT_THROW(stan::math::check_simplex("checkSimplex", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[0] = 0.5;
    y_i[1] = nan;
  }
  EXPECT_THROW(stan::math::check_simplex("checkSimplex", "y", y),
               std::domain_error);

  for (auto& y_i : y) {
    y_i[0] = nan;
  }
  EXPECT_THROW(stan::math::check_simplex("checkSimplex", "y", y),
               std::domain_error);
}
