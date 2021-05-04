#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(MathMatrixCL, optional_broadcast_errors) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
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

  // extra parenthesis needed for passing type name with comma in a macro
  EXPECT_NO_THROW((stan::math::optional_broadcast<false, false>(m_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, false>(m_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<false, true>(m_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, true>(m_cl)));

  EXPECT_NO_THROW((stan::math::optional_broadcast<false, false>(row_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, false>(row_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<false, true>(row_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, true>(row_cl)));

  EXPECT_NO_THROW((stan::math::optional_broadcast<false, false>(col_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, false>(col_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<false, true>(col_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, true>(col_cl)));

  EXPECT_NO_THROW((stan::math::optional_broadcast<false, false>(scal_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, false>(scal_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<false, true>(scal_cl)));
  EXPECT_NO_THROW((stan::math::optional_broadcast<true, true>(scal_cl)));

  EXPECT_THROW((stan::math::optional_broadcast<true, false>(col_cl) + m_cl),
               std::invalid_argument);
  EXPECT_THROW((stan::math::optional_broadcast<false, true>(row_cl) + m_cl),
               std::invalid_argument);
  EXPECT_THROW((stan::math::optional_broadcast<true, false>(scal_cl) + m_cl),
               std::invalid_argument);
}

TEST(MathMatrixCL, optional_broadcast_colwise_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(1, 3);
  m2 << 10, 100, 1000;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::colwise_optional_broadcast(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m1 + m2.replicate<3, 1>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);

  auto tmp2 = m1_cl + stan::math::colwise_optional_broadcast(m1_cl);
  matrix_cl<double> res2_cl = tmp2;
  MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);
  MatrixXd correct2 = m1 + m1;
  EXPECT_MATRIX_NEAR(res2, correct2, 1e-9);
}

TEST(MathMatrixCL, optional_broadcast_rowwise_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 1);
  m2 << 10, 100, 1000;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::rowwise_optional_broadcast(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m1 + m2.replicate<1, 3>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);

  auto tmp2 = m1_cl + stan::math::rowwise_optional_broadcast(m1_cl);
  matrix_cl<double> res2_cl = tmp2;
  MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);
  MatrixXd correct2 = m1 + m1;
  EXPECT_MATRIX_NEAR(res2, correct2, 1e-9);
}

TEST(MathMatrixCL, optional_broadcast_both_test) {
  using Eigen::MatrixXd;
  using Eigen::MatrixXi;
  using stan::math::matrix_cl;

  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(1, 1);
  m2 << 10;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + stan::math::optional_broadcast<true, true>(m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);
  MatrixXd correct = m1 + m2.replicate<3, 3>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);

  auto tmp2 = m1_cl + stan::math::optional_broadcast<true, true>(m1_cl);
  matrix_cl<double> res2_cl = tmp2;
  MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);
  MatrixXd correct2 = m1 + m1;
  EXPECT_MATRIX_NEAR(res2, correct2, 1e-9);
}

#endif
