#include <stan/math/mix.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradMixMatrixAddons, matrix_fv) {
  using Eigen::MatrixXd;
  using stan::math::matrix_fv;

  MatrixXd vals = MatrixXd::Random(100, 100);
  MatrixXd derivs = MatrixXd::Random(100, 100);

  matrix_fv mat_in = vals;

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      mat_in(i, j).d_ = derivs(i, j);
    }
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, mat_in.val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), mat_in.val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, mat_in.d().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), mat_in.d().val().array().exp());

  EXPECT_EQ(mat_in.val().val().rows(), vals.rows());
  EXPECT_EQ(mat_in.val().val().cols(), vals.cols());

  EXPECT_EQ(mat_in.d().val().rows(), derivs.rows());
  EXPECT_EQ(mat_in.d().val().cols(), derivs.cols());
}

TEST(AgradMixMatrixAddons, matrix_ffv) {
  using Eigen::MatrixXd;
  using stan::math::matrix_ffv;

  MatrixXd vals = MatrixXd::Random(100, 100);
  MatrixXd derivs = MatrixXd::Random(100, 100);

  matrix_ffv mat_in = vals;

  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 100; j++) {
      mat_in(i, j).d_ = derivs(i, j);
    }
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, mat_in.val().val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(),
                         mat_in.val().val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, mat_in.d().val().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(),
                         mat_in.d().val().val().array().exp());

  EXPECT_EQ(mat_in.val().val().val().rows(), vals.rows());
  EXPECT_EQ(mat_in.val().val().val().cols(), vals.cols());

  EXPECT_EQ(mat_in.d().val().val().rows(), derivs.rows());
  EXPECT_EQ(mat_in.d().val().val().cols(), derivs.cols());
}

TEST(AgradMixMatrixAddons, vector_fv) {
  using Eigen::VectorXd;
  using stan::math::vector_fv;

  VectorXd vals = VectorXd::Random(100);
  VectorXd derivs = VectorXd::Random(100);

  vector_fv vec_in = vals;

  for (int i = 0; i < 100; i++) {
    vec_in(i).d_ = derivs(i);
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, vec_in.val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), vec_in.val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, vec_in.d().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), vec_in.d().val().array().exp());

  EXPECT_EQ(vec_in.val().val().rows(), vals.rows());
  EXPECT_EQ(vec_in.val().val().cols(), vals.cols());

  EXPECT_EQ(vec_in.d().val().rows(), derivs.rows());
  EXPECT_EQ(vec_in.d().val().cols(), derivs.cols());
}

TEST(AgradMixMatrixAddons, vector_ffv) {
  using Eigen::VectorXd;
  using stan::math::vector_ffv;

  VectorXd vals = VectorXd::Random(100);
  VectorXd derivs = VectorXd::Random(100);

  vector_ffv vec_in = vals;

  for (int i = 0; i < 100; i++) {
    vec_in(i).d_ = derivs(i);
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, vec_in.val().val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(),
                         vec_in.val().val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, vec_in.d().val().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(),
                         vec_in.d().val().val().array().exp());

  EXPECT_EQ(vec_in.val().val().val().rows(), vals.rows());
  EXPECT_EQ(vec_in.val().val().val().cols(), vals.cols());

  EXPECT_EQ(vec_in.d().val().val().rows(), derivs.rows());
  EXPECT_EQ(vec_in.d().val().val().cols(), derivs.cols());
}

TEST(AgradMixMatrixAddons, row_vector_fv) {
  using Eigen::RowVectorXd;
  using stan::math::row_vector_fv;

  RowVectorXd vals = RowVectorXd::Random(100);
  RowVectorXd derivs = RowVectorXd::Random(100);

  row_vector_fv row_vec_in = vals;

  for (int i = 0; i < 100; i++) {
    row_vec_in(i).d_ = derivs(i);
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, row_vec_in.val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(),
                         row_vec_in.val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, row_vec_in.d().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(),
                         row_vec_in.d().val().array().exp());

  EXPECT_EQ(row_vec_in.val().val().rows(), vals.rows());
  EXPECT_EQ(row_vec_in.val().val().cols(), vals.cols());

  EXPECT_EQ(row_vec_in.d().val().rows(), derivs.rows());
  EXPECT_EQ(row_vec_in.d().val().cols(), derivs.cols());
}

TEST(AgradMixMatrixAddons, row_vector_ffv) {
  using Eigen::RowVectorXd;
  using stan::math::row_vector_ffv;

  RowVectorXd vals = RowVectorXd::Random(100);
  RowVectorXd derivs = RowVectorXd::Random(100);

  row_vector_ffv row_vec_in = vals;

  for (int i = 0; i < 100; i++) {
    row_vec_in(i).d_ = derivs(i);
  }

  EXPECT_MATRIX_FLOAT_EQ(vals, row_vec_in.val().val().val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(),
                         row_vec_in.val().val().val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, row_vec_in.d().val().val());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(),
                         row_vec_in.d().val().val().array().exp());

  EXPECT_EQ(row_vec_in.val().val().val().rows(), vals.rows());
  EXPECT_EQ(row_vec_in.val().val().val().cols(), vals.cols());

  EXPECT_EQ(row_vec_in.d().val().val().rows(), derivs.rows());
  EXPECT_EQ(row_vec_in.d().val().val().cols(), derivs.cols());
}
