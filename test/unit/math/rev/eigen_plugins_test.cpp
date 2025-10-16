#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrixAddons, var_matrix) {
  using Eigen::MatrixXd;
  using stan::math::matrix_v;
  using stan::math::matrix_vi;

  MatrixXd vals = MatrixXd::Random(100, 100);
  MatrixXd derivs = MatrixXd::Random(100, 100);

  matrix_v mat_in = vals;

  for (int i = 0; i < 100; i++)
    for (int j = 0; j < 100; j++)
      mat_in(i, j).vi_->adj_ = derivs(i, j);

  EXPECT_MATRIX_FLOAT_EQ(vals, mat_in.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.val(), mat_in.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), mat_in.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, mat_in.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), mat_in.adj().array().exp());

  EXPECT_EQ(mat_in.val().rows(), vals.rows());
  EXPECT_EQ(mat_in.val().cols(), vals.cols());

  EXPECT_EQ(mat_in.adj().rows(), derivs.rows());
  EXPECT_EQ(mat_in.adj().cols(), derivs.cols());

  const matrix_v const_mat_in = matrix_v::Random(100, 100);

  MatrixXd tri_out = const_mat_in.val().triangularView<Eigen::Upper>().solve(
      const_mat_in.adj().transpose());

  matrix_vi mat_vi = mat_in.vi();

  EXPECT_MATRIX_FLOAT_EQ(vals, mat_vi.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), mat_vi.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, mat_vi.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), mat_vi.adj().array().exp());

  EXPECT_EQ(mat_vi.val().rows(), vals.rows());
  EXPECT_EQ(mat_vi.val().cols(), vals.cols());

  EXPECT_EQ(mat_vi.adj().rows(), derivs.rows());
  EXPECT_EQ(mat_vi.adj().cols(), derivs.cols());
}

TEST(AgradRevMatrixAddons, var_vector) {
  using Eigen::VectorXd;
  using stan::math::vector_v;
  using stan::math::vector_vi;

  VectorXd vals = VectorXd::Random(100);
  VectorXd derivs = VectorXd::Random(100);

  vector_v vec_in = vals;

  for (int i = 0; i < 100; i++)
    vec_in(i).vi_->adj_ = derivs(i);

  EXPECT_MATRIX_FLOAT_EQ(vals, vec_in.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), vec_in.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, vec_in.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), vec_in.adj().array().exp());

  EXPECT_EQ(vec_in.val().rows(), vals.rows());
  EXPECT_EQ(vec_in.val().cols(), vals.cols());

  EXPECT_EQ(vec_in.adj().rows(), derivs.rows());
  EXPECT_EQ(vec_in.adj().cols(), derivs.cols());

  vector_vi vec_vi = vec_in.vi();

  EXPECT_MATRIX_FLOAT_EQ(vals, vec_vi.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), vec_vi.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, vec_vi.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), vec_vi.adj().array().exp());

  EXPECT_EQ(vec_vi.val().rows(), vals.rows());
  EXPECT_EQ(vec_vi.val().cols(), vals.cols());

  EXPECT_EQ(vec_vi.adj().rows(), derivs.rows());
  EXPECT_EQ(vec_vi.adj().cols(), derivs.cols());
}

TEST(AgradRevMatrixAddons, var_row_vector) {
  using Eigen::RowVectorXd;
  using stan::math::row_vector_v;
  using stan::math::row_vector_vi;

  RowVectorXd vals = RowVectorXd::Random(100);
  RowVectorXd derivs = RowVectorXd::Random(100);

  row_vector_v row_vec_in = vals;

  for (int i = 0; i < 100; i++)
    row_vec_in(i).vi_->adj_ = derivs(i);

  EXPECT_MATRIX_FLOAT_EQ(vals, row_vec_in.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), row_vec_in.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, row_vec_in.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), row_vec_in.adj().array().exp());

  EXPECT_EQ(row_vec_in.val().rows(), vals.rows());
  EXPECT_EQ(row_vec_in.val().cols(), vals.cols());

  EXPECT_EQ(row_vec_in.adj().rows(), derivs.rows());
  EXPECT_EQ(row_vec_in.adj().cols(), derivs.cols());

  row_vector_vi row_vec_vi = row_vec_in.vi();

  EXPECT_MATRIX_FLOAT_EQ(vals, row_vec_vi.val());
  EXPECT_MATRIX_FLOAT_EQ(vals.array().exp(), row_vec_vi.val().array().exp());

  EXPECT_MATRIX_FLOAT_EQ(derivs, row_vec_vi.adj());
  EXPECT_MATRIX_FLOAT_EQ(derivs.array().exp(), row_vec_vi.adj().array().exp());

  EXPECT_EQ(row_vec_vi.val().rows(), vals.rows());
  EXPECT_EQ(row_vec_vi.val().cols(), vals.cols());

  EXPECT_EQ(row_vec_vi.adj().rows(), derivs.rows());
  EXPECT_EQ(row_vec_vi.adj().cols(), derivs.cols());
}
