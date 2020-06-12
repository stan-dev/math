#include <stan/math/fwd.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrix, read_fvar_mat_fd) {
  using Eigen::MatrixXd;
  using stan::math::matrix_fd;
  using stan::math::read_fvar;

  matrix_fd matrix_fvar(100, 100);
  MatrixXd matrix_val(100, 100), matrix_deriv(100, 100);
  matrix_fvar.val() = MatrixXd::Random(100, 100);
  matrix_fvar.d() = MatrixXd::Random(100, 100);
  read_fvar(matrix_fvar, matrix_val, matrix_deriv);

  expect_matrix_eq(matrix_fvar.val(), matrix_val);
  expect_matrix_eq(matrix_fvar.d(), matrix_deriv);
}

TEST(AgradFwdMatrix, read_fvar_vec_fd) {
  using Eigen::VectorXd;
  using stan::math::read_fvar;
  using stan::math::vector_fd;

  vector_fd vector_fvar(100);
  VectorXd vector_val(100), vector_deriv(100);
  vector_fvar.val() = VectorXd::Random(100);
  vector_fvar.d() = VectorXd::Random(100);
  read_fvar(vector_fvar, vector_val, vector_deriv);

  expect_matrix_eq(vector_fvar.val(), vector_val);
  expect_matrix_eq(vector_fvar.d(), vector_deriv);
}

TEST(AgradFwdMatrix, read_fvar_rowvec_fd) {
  using Eigen::RowVectorXd;
  using stan::math::read_fvar;
  using stan::math::row_vector_fd;

  row_vector_fd row_vector_fvar(100);
  RowVectorXd row_vector_val(100), row_vector_deriv(100);
  row_vector_fvar.val() = RowVectorXd::Random(100);
  row_vector_fvar.d() = RowVectorXd::Random(100);
  read_fvar(row_vector_fvar, row_vector_val, row_vector_deriv);

  expect_matrix_eq(row_vector_fvar.val(), row_vector_val);
  expect_matrix_eq(row_vector_fvar.d(), row_vector_deriv);
}

TEST(AgradFwdMatrix, read_fvar_expr_fd) {
  using Eigen::MatrixXd;
  using stan::math::matrix_fd;
  using stan::math::read_fvar;

  matrix_fd matrix_fvar(100, 100);
  Eigen::VectorXd matrix_diag_val(100), matrix_diag_deriv(100);
  matrix_fvar.val() = MatrixXd::Random(100, 100);
  matrix_fvar.d() = MatrixXd::Random(100, 100);
  read_fvar(matrix_fvar.diagonal(), matrix_diag_val, matrix_diag_deriv);

  expect_matrix_eq(matrix_fvar.diagonal().val(), matrix_diag_val);
  expect_matrix_eq(matrix_fvar.diagonal().d(), matrix_diag_deriv);
}

TEST(AgradFwdMatrix, read_fvar_mat_ffd) {
  using Eigen::MatrixXd;
  using stan::math::matrix_fd;
  using stan::math::matrix_ffd;
  using stan::math::read_fvar;

  matrix_ffd matrix_fvar_fvar(100, 100);
  matrix_fd matrix_val(100, 100), matrix_deriv(100, 100);
  matrix_fvar_fvar.val().val() = MatrixXd::Random(100, 100);
  matrix_fvar_fvar.d().val() = MatrixXd::Random(100, 100);
  read_fvar(matrix_fvar_fvar, matrix_val, matrix_deriv);

  expect_matrix_eq(matrix_fvar_fvar.val().val(), matrix_val.val());
  expect_matrix_eq(matrix_fvar_fvar.d().val(), matrix_deriv.val());
}

TEST(AgradFwdMatrix, read_fvar_vec_ffd) {
  using Eigen::VectorXd;
  using stan::math::read_fvar;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  vector_ffd vector_fvar_fvar(100);
  vector_fd vector_val(100), vector_deriv(100);
  vector_fvar_fvar.val().val() = VectorXd::Random(100);
  vector_fvar_fvar.d().val() = VectorXd::Random(100);
  read_fvar(vector_fvar_fvar, vector_val, vector_deriv);

  expect_matrix_eq(vector_fvar_fvar.val().val(), vector_val.val());
  expect_matrix_eq(vector_fvar_fvar.d().val(), vector_deriv.val());
}

TEST(AgradFwdMatrix, read_fvar_rowvec_ffd) {
  using Eigen::RowVectorXd;
  using stan::math::read_fvar;
  using stan::math::row_vector_fd;
  using stan::math::row_vector_ffd;

  row_vector_ffd row_vector_fvar_fvar(100);
  row_vector_fd row_vector_val(100), row_vector_deriv(100);
  row_vector_fvar_fvar.val().val() = RowVectorXd::Random(100);
  row_vector_fvar_fvar.d().val() = RowVectorXd::Random(100);
  read_fvar(row_vector_fvar_fvar, row_vector_val, row_vector_deriv);

  expect_matrix_eq(row_vector_fvar_fvar.val().val(), row_vector_val.val());
  expect_matrix_eq(row_vector_fvar_fvar.d().val(), row_vector_deriv.val());
}

TEST(AgradFwdMatrix, read_fvar_expr_ffd) {
  using Eigen::MatrixXd;
  using stan::math::matrix_ffd;
  using stan::math::read_fvar;
  using stan::math::vector_fd;

  matrix_ffd matrix_fvar_fvar(100, 100);
  vector_fd matrix_diag_val(100), matrix_diag_deriv(100);
  matrix_fvar_fvar.val().val() = MatrixXd::Random(100, 100);
  matrix_fvar_fvar.d().val() = MatrixXd::Random(100, 100);
  read_fvar(matrix_fvar_fvar.diagonal(), matrix_diag_val, matrix_diag_deriv);

  expect_matrix_eq(matrix_fvar_fvar.diagonal().val().val(),
                   matrix_diag_val.val());
  expect_matrix_eq(matrix_fvar_fvar.diagonal().d().val(),
                   matrix_diag_deriv.val());
}
