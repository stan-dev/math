#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>

TEST(ErrorHandlingMat, CheckLess_Matrix) {
  using stan::math::check_less;
  const char* function = "check_less";
  double x = 0;
  double high = 0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_mat;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> high_mat;
  x_mat.resize(3, 3);
  high_mat.resize(3, 3);
  std::vector<double> x_scalar_vec{x, x, x};
  std::vector<double> high_scalar_vec{high, high, high};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> x_vec{
      x_mat, x_mat, x_mat};
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> high_vec{
      high_mat, high_mat, high_mat};
  // x_vec, high
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i] << -5, 0, 5, -5, 0, 5, -5, 0, 5;
    high_scalar_vec[i] = 10;
  }
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_scalar_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    high_scalar_vec[i] = std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_scalar_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, 5;
    high_scalar_vec[i] = 5;
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_scalar_vec),
               std::domain_error);

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, std::numeric_limits<double>::infinity();
    high_scalar_vec[i] = 5;
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_scalar_vec),
               std::domain_error);

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, std::numeric_limits<double>::infinity();
    high_scalar_vec[i] = std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_scalar_vec),
               std::domain_error);

  // x_vec, high_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i] << -5, 0, 5, -5, 0, 5, -5, 0, 5;
    high_vec[i] << 0, 5, 10, 0, 5, 10, 0, 5, 10;
  }
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, 5;
    high_vec[i].col(i) << std::numeric_limits<double>::infinity(), 10, 10;
  }
  EXPECT_NO_THROW(check_less(function, "x", x_vec, high_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, 5;
    high_vec[i].col(i) << 10, 10, 5;
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec), std::domain_error);

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, std::numeric_limits<double>::infinity();
    high_vec[i].col(i) << 10, 10, 10;
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec), std::domain_error);

  for (int i = 0; i < x_vec.size(); ++i) {
    x_vec[i].col(i) << -5, 0, std::numeric_limits<double>::infinity();
    high_vec[i].col(i) << 10, 10, std::numeric_limits<double>::infinity();
  }
  EXPECT_THROW(check_less(function, "x", x_vec, high_vec), std::domain_error);

  // x, high_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = -100;
    high_vec[i].col(i) << 0, 5, 10;
  }
  EXPECT_NO_THROW(check_less(function, "x", x_scalar_vec, high_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 10;
    high_vec[i] << 100, 200, std::numeric_limits<double>::infinity(), 100, 200,
        std::numeric_limits<double>::infinity(), 100, 200,
        std::numeric_limits<double>::infinity();
  }
  EXPECT_NO_THROW(check_less(function, "x", x_scalar_vec, high_vec));

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 5;
    high_vec[i] << 100, 200, 5, 100, 200, 5, 100, 200, 5;
  }
  EXPECT_THROW(check_less(function, "x", x_scalar_vec, high_vec),
               std::domain_error);

  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = std::numeric_limits<double>::infinity();
    high_vec[i] << 10, 20, 30, 10, 20, 30, 10, 20, 30;
  }
  EXPECT_THROW(check_less(function, "x", x_scalar_vec, high_vec),
               std::domain_error);

  constexpr auto inf_val = std::numeric_limits<double>::infinity();
  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = inf_val;
    high_vec[i] << inf_val, inf_val, inf_val, inf_val, inf_val, inf_val,
        inf_val, inf_val, inf_val;
  }
  EXPECT_THROW(check_less(function, "x", x_scalar_vec, high_vec),
               std::domain_error);
}

TEST(ErrorHandlingMat, CheckLess_Matrix_one_indexed_message) {
  using stan::math::check_less;
  const char* function = "check_less";
  double x = 0;
  double high = 0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_mat;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> high_mat;
  x_mat.resize(3, 3);
  high_mat.resize(3, 3);
  std::vector<double> x_scalar_vec{x, x, x};
  std::vector<double> high_scalar_vec{high, high, high};
  using std_vec_mat
      = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>;
  std_vec_mat x_vec1{x_mat, x_mat, x_mat};
  std_vec_mat high_vec1{high_mat, high_mat, high_mat};
  // x_vec, high
  for (int i = 0; i < x_vec1.size(); ++i) {
    x_vec1[i] << -5, 0, 5, -5, 0, 5, -5, 0, 5;
    high_scalar_vec[i] = 5;
  }
  std::vector<std_vec_mat> x_vec{x_vec1, x_vec1, x_vec1};
  std::vector<std_vec_mat> high_vec{high_vec1, high_vec1, high_vec1};
  std::string message;
  try {
    check_less(function, "x", x_vec, high_scalar_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("[1, 1][1, 3]")) << message;

  // x_vec, high_vec
  for (auto& x_vec_i : x_vec) {
    for (auto& x_vec1_i : x_vec_i) {
      x_vec1_i << -5, 5, 0, -5, 5, 0, -5, 5, 0;
    }
  }

  for (auto& high_vec_i : high_vec) {
    for (auto& high_vec1_i : high_vec_i) {
      high_vec1_i << -5, 5, 0, -5, 5, 0, -5, 5, 0;
    }
  }

  try {
    check_less(function, "x", x_vec, high_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }
  EXPECT_NE(std::string::npos, message.find("[1, 1][1, 1]")) << message;

  // x, high_vec
  for (int i = 0; i < x_vec.size(); ++i) {
    x_scalar_vec[i] = 30;
    high_vec[i][i].col(i) << 10, 20, 30;
  }
  try {
    check_less(function, "x", x_scalar_vec, high_vec);
    FAIL() << "should have thrown";
  } catch (std::domain_error& e) {
    message = e.what();
  } catch (...) {
    FAIL() << "threw the wrong error";
  }

  EXPECT_NE(std::string::npos, message.find("["))
      << "index provided when x has none" << std::endl
      << message;
}

TEST(ErrorHandlingMat, CheckLess_nan) {
  using stan::math::check_less;
  const char* function = "check_less";
  double nan = std::numeric_limits<double>::quiet_NaN();

  Eigen::Matrix<double, Eigen::Dynamic, 1> x_vec(3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> low_vec(3);

  // x_vec, low_vec
  x_vec << -1, 0, 1;
  low_vec << -2, -1, 0;
  EXPECT_THROW(check_less(function, "x", x_vec, nan), std::domain_error);

  for (int i = 0; i < x_vec.size(); i++) {
    x_vec << -1, 0, 1;
    x_vec(i) = nan;
    EXPECT_THROW(check_less(function, "x", x_vec, low_vec), std::domain_error);
  }

  x_vec << -1, 0, 1;
  for (int i = 0; i < low_vec.size(); i++) {
    low_vec << -1, 0, 1;
    low_vec(i) = nan;
    EXPECT_THROW(check_less(function, "x", x_vec, low_vec), std::domain_error);
  }

  for (int i = 0; i < x_vec.size(); i++) {
    x_vec << -1, 0, 1;
    low_vec << -2, -1, 0;
    x_vec(i) = nan;
    for (int j = 0; j < low_vec.size(); j++) {
      low_vec(i) = nan;
      EXPECT_THROW(check_less(function, "x", x_vec, low_vec),
                   std::domain_error);
    }
  }
}

TEST(ErrorHandlingScalar, CheckLess) {
  using stan::math::check_less;
  const char* function = "check_less";
  double x = -10.0;
  double lb = 0.0;

  EXPECT_NO_THROW(check_less(function, "x", x, lb))
      << "check_less should be true with x < lb";

  x = 1.0;
  EXPECT_THROW(check_less(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x > lb";

  x = lb;
  EXPECT_THROW(check_less(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x == lb";

  x = -std::numeric_limits<double>::infinity();
  EXPECT_NO_THROW(check_less(function, "x", x, lb))
      << "check_less should be true with x == -Inf and lb = 0.0";

  x = -10.0;
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x == -10.0 and lb == -Inf";

  x = -std::numeric_limits<double>::infinity();
  lb = -std::numeric_limits<double>::infinity();
  EXPECT_THROW(check_less(function, "x", x, lb), std::domain_error)
      << "check_less should throw an exception with x == -Inf and lb == -Inf";
}

TEST(ErrorHandlingScalar, CheckLess_nan) {
  using stan::math::check_less;
  const char* function = "check_less";
  double x = 10.0;
  double lb = 0.0;
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(check_less(function, "x", nan, lb), std::domain_error);
  EXPECT_THROW(check_less(function, "x", x, nan), std::domain_error);
  EXPECT_THROW(check_less(function, "x", nan, nan), std::domain_error);
}
