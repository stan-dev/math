#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include <vector>

using stan::math::var;

TEST(AgradRev_count_vars, int_arg) {
  int arg = 5;

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, double_arg) {
  double arg = 5.0;

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_int_arg) {
  std::vector<int> arg(5, 10);

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_double_arg) {
  std::vector<double> arg(5, 10.0);

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_vector_arg) {
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_row_vector_arg) {
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_matrix_arg) {
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_std_vector_double_arg) {
  std::vector<std::vector<double>> arg(5, std::vector<double>(5, 10.0));

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_vector_arg) {
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(5));

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_row_vector_arg) {
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::VectorXd::Ones(5));

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_matrix_arg) {
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(5, 3));

  EXPECT_EQ(0, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, var_arg) {
  var arg = 5.0;

  EXPECT_EQ(1, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_var_arg) {
  std::vector<var> arg(5);

  EXPECT_EQ(5, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_vector_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg(5);

  EXPECT_EQ(5, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_row_vector_var_arg) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg(5);

  EXPECT_EQ(5, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, eigen_matrix_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg(5, 5);

  EXPECT_EQ(25, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_std_vector_var_arg) {
  std::vector<std::vector<var>> arg(5, std::vector<var>(5));

  EXPECT_EQ(25, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_vector_var_arg) {
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> arg(
      2, Eigen::Matrix<var, Eigen::Dynamic, 1>(5));

  EXPECT_EQ(10, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_row_vector_var_arg) {
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> arg(
      2, Eigen::Matrix<var, 1, Eigen::Dynamic>(5));

  EXPECT_EQ(10, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, std_vector_eigen_matrix_var_arg) {
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arg(
      2, Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>(5, 3));

  EXPECT_EQ(30, stan::math::count_vars(arg));
}

TEST(AgradRev_count_vars, zero_args) { EXPECT_EQ(0, stan::math::count_vars()); }

TEST(AgradRev_count_vars, sum) {
  int arg1 = 1;
  double arg2 = 1.0;
  std::vector<int> arg3(5, 1);
  std::vector<double> arg4(5, 1.0);
  Eigen::VectorXd arg5 = Eigen::VectorXd::Ones(5);
  Eigen::RowVectorXd arg6 = Eigen::RowVectorXd::Ones(5);
  Eigen::MatrixXd arg7 = Eigen::MatrixXd::Ones(5, 5);
  std::vector<std::vector<double>> arg8(2, arg4);
  std::vector<Eigen::VectorXd> arg9(2, arg5);

  var arg10 = 1.0;
  std::vector<var> arg11(5, 1.0);
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg12(3);
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg13(4);
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg14(5, 3);
  std::vector<std::vector<var>> arg15(2, arg11);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> arg16(2, arg12);
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> arg17(2, arg13);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arg18(2,
                                                                        arg14);

  EXPECT_EQ(
      1 + 5 + 3 + 4 + 15 + 2 * 5 + 2 * 3 + 2 * 4 + 30,
      count_vars(arg1, arg18, arg17, arg2, arg16, arg3, arg15, arg4, arg14,
                 arg5, arg13, arg12, arg6, arg11, arg7, arg10, arg8, arg9));
}
