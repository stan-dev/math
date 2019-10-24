#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixAddons, fvar_double_matrix) {
  using Eigen::MatrixXd;
  using stan::math::matrix_fd;

  MatrixXd vals = MatrixXd::Random(100, 100);
  MatrixXd derivs = MatrixXd::Random(100, 100);

  matrix_fd mat_in(100, 100);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      mat_in(i, j).val_ = vals(i, j);
      mat_in(i, j).d_ = derivs(i, j);
    }
  }

  expect_matrix_eq(vals, mat_in.val());
  expect_matrix_eq(vals.array().exp(), mat_in.val().array().exp());

  expect_matrix_eq(derivs, mat_in.d());
  expect_matrix_eq(derivs.array().exp(), mat_in.d().array().exp());

  EXPECT_EQ(mat_in.val().rows(), vals.rows());
  EXPECT_EQ(mat_in.val().cols(), vals.cols());

  EXPECT_EQ(mat_in.d().rows(), derivs.rows());
  EXPECT_EQ(mat_in.d().cols(), derivs.cols());
}

TEST(AgradFwdMatrixAddons, fvarfvar_double_matrix) {
  using Eigen::MatrixXd;
  using stan::math::matrix_ffd;

  MatrixXd vals = MatrixXd::Random(100, 100);
  MatrixXd derivs = MatrixXd::Random(100, 100);

  matrix_ffd mat_in(100, 100);

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      mat_in(i, j).val_.val_ = vals(i, j);
      mat_in(i, j).d_.val_ = derivs(i, j);
    }
  }

  expect_matrix_eq(vals, mat_in.val().val());
  expect_matrix_eq(vals.array().exp(), mat_in.val().val().array().exp());

  expect_matrix_eq(derivs, mat_in.d().val());
  expect_matrix_eq(derivs.array().exp(), mat_in.d().val().array().exp());

  EXPECT_EQ(mat_in.val().rows(), vals.rows());
  EXPECT_EQ(mat_in.val().cols(), vals.cols());

  EXPECT_EQ(mat_in.d().rows(), derivs.rows());
  EXPECT_EQ(mat_in.d().cols(), derivs.cols());
}

TEST(AgradFwdMatrixAddons, fvar_double_vector) {
  using Eigen::VectorXd;
  using stan::math::vector_fd;

  VectorXd vals = VectorXd::Random(100);
  VectorXd derivs = VectorXd::Random(100);

  vector_fd vec_in(100);

  for (int i = 0; i < 100; i++) {
    vec_in(i).val_ = vals(i);
    vec_in(i).d_ = derivs(i);
  }

  expect_matrix_eq(vals, vec_in.val());
  expect_matrix_eq(vals.array().exp(), vec_in.val().array().exp());

  expect_matrix_eq(derivs, vec_in.d());
  expect_matrix_eq(derivs.array().exp(), vec_in.d().array().exp());

  EXPECT_EQ(vec_in.val().rows(), vals.rows());
  EXPECT_EQ(vec_in.val().cols(), vals.cols());

  EXPECT_EQ(vec_in.d().rows(), derivs.rows());
  EXPECT_EQ(vec_in.d().cols(), derivs.cols());
}

TEST(AgradFwdMatrixAddons, fvarfvar_double_vector) {
  using Eigen::VectorXd;
  using stan::math::vector_ffd;

  VectorXd vals = VectorXd::Random(100);
  VectorXd derivs = VectorXd::Random(100);

  vector_ffd vec_in(100);

  for (int i = 0; i < 100; i++) {
    vec_in(i).val_.val_ = vals(i);
    vec_in(i).d_.val_ = derivs(i);
  }

  expect_matrix_eq(vals, vec_in.val().val());
  expect_matrix_eq(vals.array().exp(), vec_in.val().val().array().exp());

  expect_matrix_eq(derivs, vec_in.d().val());
  expect_matrix_eq(derivs.array().exp(), vec_in.d().val().array().exp());

  EXPECT_EQ(vec_in.val().rows(), vals.rows());
  EXPECT_EQ(vec_in.val().cols(), vals.cols());

  EXPECT_EQ(vec_in.d().rows(), derivs.rows());
  EXPECT_EQ(vec_in.d().cols(), derivs.cols());
}

TEST(AgradFwdMatrixAddons, fvar_double_rowvector) {
  using Eigen::RowVectorXd;
  using stan::math::row_vector_fd;

  RowVectorXd vals = RowVectorXd::Random(100);
  RowVectorXd derivs = RowVectorXd::Random(100);

  row_vector_fd row_vec_in(100);

  for (int i = 0; i < 100; i++) {
    row_vec_in(i).val_ = vals(i);
    row_vec_in(i).d_ = derivs(i);
  }

  expect_matrix_eq(vals, row_vec_in.val());
  expect_matrix_eq(vals.array().exp(), row_vec_in.val().array().exp());

  expect_matrix_eq(derivs, row_vec_in.d());
  expect_matrix_eq(derivs.array().exp(), row_vec_in.d().array().exp());

  EXPECT_EQ(row_vec_in.val().rows(), vals.rows());
  EXPECT_EQ(row_vec_in.val().cols(), vals.cols());

  EXPECT_EQ(row_vec_in.d().rows(), derivs.rows());
  EXPECT_EQ(row_vec_in.d().cols(), derivs.cols());
}

TEST(AgradFwdMatrixAddons, fvarfvar_double_rowvector) {
  using Eigen::RowVectorXd;
  using stan::math::row_vector_ffd;

  RowVectorXd vals = RowVectorXd::Random(100);
  RowVectorXd derivs = RowVectorXd::Random(100);

  row_vector_ffd row_vec_in(100);

  for (int i = 0; i < 100; i++) {
    row_vec_in(i).val_.val_ = vals(i);
    row_vec_in(i).d_.val_ = derivs(i);
  }

  expect_matrix_eq(vals, row_vec_in.val().val());
  expect_matrix_eq(vals.array().exp(), row_vec_in.val().val().array().exp());

  expect_matrix_eq(derivs, row_vec_in.d().val());
  expect_matrix_eq(derivs.array().exp(), row_vec_in.d().val().array().exp());

  EXPECT_EQ(row_vec_in.val().rows(), vals.rows());
  EXPECT_EQ(row_vec_in.val().cols(), vals.cols());

  EXPECT_EQ(row_vec_in.d().rows(), derivs.rows());
  EXPECT_EQ(row_vec_in.d().cols(), derivs.cols());
}
