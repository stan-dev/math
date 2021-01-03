#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(MathMatrixCL, broadcast_errors) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;

  MatrixXd m(3, 3);
  m << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd row(1, 3);
  row << 10, 100, 1000;
  MatrixXd col(3, 1);
  col << 10, 100, 1000;
  MatrixXd scal(1, 1);
  scal << 10;

  matrix_cl<double> m_cl(m);
  matrix_cl<double> row_cl(row);
  matrix_cl<double> col_cl(col);
  matrix_cl<double> scal_cl(scal);

  // extra parenthesis is needed for passing type name with comma in a macro
  EXPECT_NO_THROW((stan::math::broadcast<false, false>(m_cl)));
  EXPECT_THROW((stan::math::broadcast<true, false>(m_cl)),
               std::invalid_argument);
  EXPECT_THROW((stan::math::broadcast<false, true>(m_cl)),
               std::invalid_argument);
  EXPECT_THROW((stan::math::colwise_broadcast(m_cl)), std::invalid_argument);
  EXPECT_THROW((stan::math::rowwise_broadcast(m_cl)), std::invalid_argument);
  EXPECT_THROW((stan::math::broadcast<true, true>(m_cl)),
               std::invalid_argument);

  EXPECT_NO_THROW((stan::math::broadcast<false, false>(row_cl)));
  EXPECT_NO_THROW((stan::math::broadcast<true, false>(row_cl)));
  EXPECT_THROW((stan::math::broadcast<false, true>(row_cl)),
               std::invalid_argument);
  EXPECT_NO_THROW((stan::math::colwise_broadcast(row_cl)));
  EXPECT_THROW((stan::math::rowwise_broadcast(row_cl)), std::invalid_argument);
  EXPECT_THROW((stan::math::broadcast<true, true>(row_cl)),
               std::invalid_argument);

  EXPECT_NO_THROW((stan::math::broadcast<false, false>(col_cl)));
  EXPECT_THROW((stan::math::broadcast<true, false>(col_cl)),
               std::invalid_argument);
  EXPECT_NO_THROW((stan::math::broadcast<false, true>(col_cl)));
  EXPECT_THROW((stan::math::colwise_broadcast(col_cl)), std::invalid_argument);
  EXPECT_NO_THROW((stan::math::rowwise_broadcast(col_cl)));
  EXPECT_THROW((stan::math::broadcast<true, true>(col_cl)),
               std::invalid_argument);

  EXPECT_NO_THROW((stan::math::broadcast<false, false>(scal_cl)));
  EXPECT_NO_THROW((stan::math::broadcast<true, false>(scal_cl)));
  EXPECT_NO_THROW((stan::math::broadcast<false, true>(scal_cl)));
  EXPECT_NO_THROW((stan::math::colwise_broadcast(scal_cl)));
  EXPECT_NO_THROW((stan::math::rowwise_broadcast(scal_cl)));
  EXPECT_NO_THROW((stan::math::broadcast<true, true>(scal_cl)));

  EXPECT_NO_THROW((stan::math::broadcast<false, false>(scal_cl).eval()));
  EXPECT_THROW((stan::math::broadcast<true, false>(scal_cl).eval()),
               std::invalid_argument);
  EXPECT_THROW((stan::math::broadcast<false, true>(scal_cl).eval()),
               std::invalid_argument);
  EXPECT_THROW((stan::math::broadcast<true, true>(scal_cl).eval()),
               std::invalid_argument);
}

TEST(MathMatrixCL, broadcast_rows_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(1, 3);
  m2 << 10, 100, 1000;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::broadcast<true, false>(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 + m2.replicate<3, 1>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, broadcast_cols_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 1);
  m2 << 10, 100, 1000;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::broadcast<false, true>(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 + m2.replicate<1, 3>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, broadcast_both_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(1, 1);
  m2 << 10;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::broadcast<true, true>(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 + m2.replicate<3, 3>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, rowwise_sum_broadcast_test) {
  using Eigen::MatrixXd;
  using stan::math::matrix_cl;
  MatrixXd m(3, 2);
  m << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6;

  matrix_cl<double> m_cl(m);

  matrix_cl<double> res_cl
      = stan::math::rowwise_broadcast(stan::math::rowwise_sum(m_cl)) + m_cl;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m.rowwise().sum().replicate(1, 2) + m;
  EXPECT_EQ(correct.rows(), res.rows());
  EXPECT_EQ(correct.cols(), res.cols());
  EXPECT_MATRIX_NEAR(correct, res, 1e-9);
}

TEST(MathMatrixCL, broadcast_view_test) {
  using stan::math::matrix_cl;

  matrix_cl<double> v_cl(2, 1);
  matrix_cl<double> m_cl(2, 2, stan::math::matrix_cl_view::Diagonal);
  matrix_cl<double> dst_block_cl(3, 3, stan::math::matrix_cl_view::Diagonal);
  matrix_cl<double> dst_cl = stan::math::broadcast<false, false>(
      block_zero_based(dst_block_cl, 0, 0, 2, 2));
  EXPECT_EQ(dst_cl.view(), stan::math::matrix_cl_view::Diagonal);
  dst_cl = stan::math::broadcast<false, true>(v_cl) + m_cl;
  EXPECT_EQ(dst_cl.view(), stan::math::matrix_cl_view::Entire);
  block_zero_based(dst_block_cl, 0, 0, 2, 2)
      = stan::math::broadcast<false, false>(
          block_zero_based(dst_block_cl, 0, 0, 2, 2));
  EXPECT_EQ(dst_block_cl.view(), stan::math::matrix_cl_view::Diagonal);
  block_zero_based(dst_block_cl, 0, 1, 2, 2)
      = stan::math::broadcast<false, true>(v_cl) + m_cl;
  EXPECT_EQ(dst_block_cl.view(), stan::math::matrix_cl_view::Upper);
  block_zero_based(dst_block_cl, 0, 1, 2, 2)
      = stan::math::broadcast<true, false>(stan::math::transpose(v_cl)) + m_cl;
  EXPECT_EQ(dst_block_cl.view(), stan::math::matrix_cl_view::Upper);
  block_zero_based(dst_block_cl, 0, 0, 2, 2)
      = stan::math::broadcast<false, true>(v_cl) + m_cl;
  EXPECT_EQ(dst_block_cl.view(), stan::math::matrix_cl_view::Entire);
}

#endif
