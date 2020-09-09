#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, read_var_mat) {
  using Eigen::MatrixXd;
  using stan::math::matrix_v;
  using stan::math::matrix_vi;
  using stan::math::read_val_adj;
  using stan::math::read_vi_adj;
  using stan::math::read_vi_val;
  using stan::math::read_vi_val_adj;

  matrix_v matrix_var(100, 100);
  matrix_vi matrix_vari(100, 100);
  MatrixXd matrix_val(100, 100), matrix_deriv(100, 100);
  matrix_var = MatrixXd::Random(100, 100);
  matrix_var.adj() = MatrixXd::Random(100, 100);

  read_vi_val_adj(matrix_var, matrix_vari, matrix_val, matrix_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.val(), matrix_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.adj(), matrix_deriv);

  matrix_val.setZero();
  matrix_deriv.setZero();
  read_val_adj(matrix_var, matrix_val, matrix_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.val(), matrix_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.adj(), matrix_deriv);

  matrix_val.setZero();
  matrix_deriv.setZero();
  read_val_adj(matrix_vari, matrix_val, matrix_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_vari.val(), matrix_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_vari.adj(), matrix_deriv);

  matrix_val.setZero();
  matrix_vi matrix_vi2(100, 100);
  read_vi_val(matrix_var, matrix_vi2, matrix_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.val(), matrix_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.val(), matrix_vi2.val());
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.adj(), matrix_vi2.adj());

  matrix_deriv.setZero();
  matrix_vi matrix_vi3(100, 100);
  read_vi_adj(matrix_var, matrix_vi3, matrix_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.adj(), matrix_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.val(), matrix_vi3.val());
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.adj(), matrix_vi3.adj());
}

TEST(AgradRevMatrix, read_var_vec) {
  using Eigen::VectorXd;
  using stan::math::read_val_adj;
  using stan::math::read_vi_adj;
  using stan::math::read_vi_val;
  using stan::math::read_vi_val_adj;
  using stan::math::vector_v;
  using stan::math::vector_vi;

  vector_v vector_var(100);
  vector_vi vector_vari(100);
  VectorXd vector_val(100), vector_deriv(100);
  vector_var = VectorXd::Random(100);
  vector_var.adj() = VectorXd::Random(100);

  read_vi_val_adj(vector_var, vector_vari, vector_val, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.adj(), vector_deriv);

  vector_val.setZero();
  vector_deriv.setZero();
  read_val_adj(vector_vari, vector_val, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(vector_vari.val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(vector_vari.adj(), vector_deriv);

  vector_val.setZero();
  vector_vi vector_vi2(100);
  read_vi_val(vector_var, vector_vi2, vector_val);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.val(), vector_vi2.val());
  EXPECT_MATRIX_FLOAT_EQ(vector_var.adj(), vector_vi2.adj());

  vector_deriv.setZero();
  vector_vi vector_vi3(100);
  read_vi_adj(vector_var, vector_vi3, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.adj(), vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(vector_var.val(), vector_vi3.val());
  EXPECT_MATRIX_FLOAT_EQ(vector_var.adj(), vector_vi3.adj());
}

TEST(AgradRevMatrix, read_var_rowvec) {
  using Eigen::RowVectorXd;
  using stan::math::read_val_adj;
  using stan::math::read_vi_adj;
  using stan::math::read_vi_val;
  using stan::math::read_vi_val_adj;
  using stan::math::row_vector_v;
  using stan::math::row_vector_vi;

  row_vector_v row_vector_var(100);
  row_vector_vi row_vector_vari(100);
  RowVectorXd row_vector_val(100), row_vector_deriv(100);
  row_vector_var = RowVectorXd::Random(100);
  row_vector_var.adj() = RowVectorXd::Random(100);

  read_vi_val_adj(row_vector_var, row_vector_vari, row_vector_val,
                  row_vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.val(), row_vector_val);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.adj(), row_vector_deriv);

  row_vector_val.setZero();
  row_vector_deriv.setZero();
  read_val_adj(row_vector_vari, row_vector_val, row_vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_vari.val(), row_vector_val);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_vari.adj(), row_vector_deriv);

  row_vector_val.setZero();
  row_vector_vi row_vector_vi2(100);
  read_vi_val(row_vector_var, row_vector_vi2, row_vector_val);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.val(), row_vector_val);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.val(), row_vector_vi2.val());
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.adj(), row_vector_vi2.adj());

  row_vector_deriv.setZero();
  row_vector_vi row_vector_vi3(100);
  read_vi_adj(row_vector_var, row_vector_vi3, row_vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.adj(), row_vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.val(), row_vector_vi3.val());
  EXPECT_MATRIX_FLOAT_EQ(row_vector_var.adj(), row_vector_vi3.adj());
}

TEST(AgradRevMatrix, read_var_expr) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using stan::math::matrix_v;
  using stan::math::matrix_vi;
  using stan::math::read_val_adj;
  using stan::math::read_vi_adj;
  using stan::math::read_vi_val;
  using stan::math::read_vi_val_adj;
  using stan::math::vector_vi;

  matrix_v matrix_var(100, 100);
  matrix_vi matrix_vari(100, 100);
  vector_vi vector_vari(100);
  VectorXd vector_val(100), vector_deriv(100);
  matrix_var = MatrixXd::Random(100, 100);
  matrix_var.adj() = MatrixXd::Random(100, 100);
  matrix_vari = matrix_var.vi();

  read_vi_val_adj(matrix_var.diagonal(), vector_vari, vector_val, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_deriv);

  vector_val.setZero();
  vector_deriv.setZero();
  read_val_adj(matrix_var.diagonal(), vector_val, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_deriv);

  vector_val.setZero();
  vector_deriv.setZero();
  read_val_adj(matrix_vari.diagonal(), vector_val, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_deriv);

  vector_val.setZero();
  vector_vi vector_vari2(100);
  read_vi_val(matrix_var.diagonal(), vector_vari2, vector_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_val);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_vari2.val());
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_vari2.adj());

  vector_deriv.setZero();
  vector_vi vector_vari3(100);
  read_vi_adj(matrix_var.diagonal(), vector_vari3, vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_deriv);
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().val(), vector_vari3.val());
  EXPECT_MATRIX_FLOAT_EQ(matrix_var.diagonal().adj(), vector_vari3.adj());
}
