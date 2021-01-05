#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, PromoteVarMatrix) {
  using stan::promote_var_matrix_t;
  using stan::math::var;
  using stan::math::var_value;
  using std::is_same;
  using var_matrix = var_value<Eigen::MatrixXd>;
  using var_vector = var_value<Eigen::VectorXd>;
  using var_row_vector = var_value<Eigen::RowVectorXd>;
  using matrix_var = Eigen::Matrix<var, -1, -1>;
  using vector_var = Eigen::Matrix<var, -1, 1>;
  using row_vector_var = Eigen::Matrix<var, 1, -1>;
  var_matrix A_vm(Eigen::MatrixXd::Zero(10, 10));
  matrix_var A_mv(Eigen::MatrixXd::Zero(10, 10));

  EXPECT_TRUE((is_same<var_value<Eigen::MatrixXd>,
                       promote_var_matrix_t<Eigen::MatrixXd, var_matrix,
                                            matrix_var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::MatrixXd>,
                       promote_var_matrix_t<Eigen::MatrixXd, var_matrix,
                                            vector_var, var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::MatrixXd>,
                       promote_var_matrix_t<Eigen::MatrixXd, var_vector,
                                            vector_var, double>>::value));

  EXPECT_TRUE((is_same<var_value<Eigen::VectorXd>,
                       promote_var_matrix_t<Eigen::VectorXd, var_matrix,
                                            matrix_var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::VectorXd>,
                       promote_var_matrix_t<Eigen::VectorXd, var_matrix,
                                            vector_var, var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::VectorXd>,
                       promote_var_matrix_t<Eigen::VectorXd, var_vector,
                                            vector_var, double>>::value));

  EXPECT_TRUE(
      (is_same<stan::math::var_value<stan::math::promote_scalar_t<
                   double, decltype(A_vm.val().block(0, 0, 2, 2))>>,
               promote_var_matrix_t<decltype(A_vm.block(0, 0, 2, 2))>>::value));

  EXPECT_TRUE(
      (is_same<matrix_var,
               promote_var_matrix_t<decltype(A_mv.block(0, 0, 2, 2))>>::value));

  EXPECT_TRUE((
      is_same<matrix_var, promote_var_matrix_t<decltype(A_mv * A_mv)>>::value));

  EXPECT_TRUE((is_same<var_value<Eigen::RowVectorXd>,
                       promote_var_matrix_t<Eigen::RowVectorXd, var_matrix,
                                            matrix_var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::RowVectorXd>,
                       promote_var_matrix_t<Eigen::RowVectorXd, var_matrix,
                                            vector_var, var>>::value));
  EXPECT_TRUE((is_same<var_value<Eigen::RowVectorXd>,
                       promote_var_matrix_t<Eigen::RowVectorXd, var_vector,
                                            row_vector_var, double>>::value));

  EXPECT_TRUE((is_same<Eigen::Matrix<var, -1, -1>,
                       promote_var_matrix_t<Eigen::MatrixXd, vector_var,
                                            vector_var, double>>::value));
  EXPECT_TRUE((is_same<Eigen::Matrix<var, 1, -1>,
                       promote_var_matrix_t<Eigen::RowVectorXd, vector_var,
                                            row_vector_var, double>>::value));
}
