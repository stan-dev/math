#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathFunctions, reverse_array) {
  using stan::math::reverse;

  std::vector<double> x0;
  EXPECT_EQ(0, reverse(x0).size());

  std::vector<double> x1(5);
  EXPECT_FLOAT_EQ(x1[0], reverse(x1)[0]);

  std::vector<double> x{7, -1.2, 4, -3};
  std::vector<double> rev{-3, 4, -1.2, 7};
  std::vector<double> y = stan::math::reverse(x);
  for (int i = 0; i < x.size(); i++) {
    EXPECT_FLOAT_EQ(y[i], rev[i]);
  }

  std::vector<std::vector<double>> cont_x{x0, x1, x};
  std::vector<std::vector<double>> cont_rev{x, x1, x0};
  std::vector<std::vector<double>> cont_y = stan::math::reverse(cont_x);
  for (int i = 0; i < cont_x.size(); i++) {
    EXPECT_EQ(cont_y[i].size(), cont_rev[i].size());
    for (int j = 0; j < cont_y[i].size(); j++) {
      EXPECT_FLOAT_EQ(cont_y[i][j], cont_rev[i][j]);
    }
  }
}

TEST(MathFunctions, reverse_vector) {
  using stan::math::reverse;

  Eigen::VectorXd x0(0);
  EXPECT_EQ(0, reverse(x0).size());

  Eigen::VectorXd x1 = Eigen::VectorXd::Constant(1, 5);
  EXPECT_MATRIX_FLOAT_EQ(x1, reverse(x1));

  Eigen::VectorXd x(4);
  x << 7, -1.2, 4, -3;
  Eigen::VectorXd rev(4);
  rev << -3, 4, -1.2, 7;
  EXPECT_MATRIX_FLOAT_EQ(rev, reverse(x));

  std::vector<Eigen::VectorXd> cont_x{x0, x1, x};
  std::vector<Eigen::VectorXd> cont_rev{x, x1, x0};
  std::vector<Eigen::VectorXd> cont_y = stan::math::reverse(cont_x);
  for (int i = 0; i < cont_x.size(); i++) {
    EXPECT_MATRIX_FLOAT_EQ(cont_y[i], cont_rev[i]);
  }
}

TEST(MathFunctions, reverse_row_vector) {
  using stan::math::reverse;

  Eigen::RowVectorXd x0(0);
  EXPECT_EQ(0, reverse(x0).size());

  Eigen::RowVectorXd x1 = Eigen::RowVectorXd::Constant(1, 5);
  EXPECT_MATRIX_FLOAT_EQ(x1, reverse(x1));

  Eigen::RowVectorXd x(4);
  x << 7, -1.2, 4, -3;
  Eigen::RowVectorXd rev(4);
  rev << -3, 4, -1.2, 7;
  EXPECT_MATRIX_FLOAT_EQ(rev, reverse(x));

  std::vector<Eigen::RowVectorXd> cont_x{x0, x1, x};
  std::vector<Eigen::RowVectorXd> cont_rev{x, x1, x0};
  std::vector<Eigen::RowVectorXd> cont_y = stan::math::reverse(cont_x);
  for (int i = 0; i < cont_x.size(); i++) {
    EXPECT_MATRIX_FLOAT_EQ(cont_y[i], cont_rev[i]);
  }
}
