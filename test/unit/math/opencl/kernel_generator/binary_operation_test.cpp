#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <string>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrixCL, addition_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + m2_cl;
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 + m2.cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, subtraction_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);

  auto tmp = m1_cl - m2_cl;
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 - m2.cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, elewise_multiplication_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);

  auto tmp = elewise_multiplication(m1_cl, m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.array() * m2.cast<double>().array();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, elewise_division_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 1, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  auto tmp = elewise_division(m1_cl, m2_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.array() / m2.cast<double>().array();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, multiple_operations) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;
  MatrixXi m3(3, 3);
  m3 << 1, 10, 1100, -40, -14, -1, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  matrix_cl<int> m3_cl(m3);
  auto tmp = m1_cl * 2. - (m2_cl + m3_cl) * 4;
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 * 2. - ((m2 + m3) * 4).cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, multiple_operations_accepts_lvalue) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;
  MatrixXi m3(3, 3);
  m3 << 1, 10, 1100, -40, -14, -1, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  matrix_cl<int> m3_cl(m3);
  auto tmp = (m2_cl + m3_cl) * 4;
  auto tmp2 = m1_cl * 2. - tmp;
  matrix_cl<double> res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 * 2. - ((m2 + m3) * 4).cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, multiplication_with_scalar_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl * 2;
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  EXPECT_MATRIX_NEAR(res, (m1.array() * 2).matrix().eval(), 1e-9);

  auto tmp2 = 2 * m1_cl;
  matrix_cl<double> res2_cl = tmp2;
  MatrixXd res2 = stan::math::from_matrix_cl(res2_cl);

  MatrixXd correct = m1.array() * 2;
  EXPECT_MATRIX_NEAR(res2, correct, 1e-9);
}

TEST(MathMatrixCL, matrix_multiplication_in_expression_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = ((m1_cl - 2) * (m2_cl + 3)) - 1;
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct
      = ((m1.array() - 2.).matrix() * (m2.array() + 3.).matrix()).array() - 1.;
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(MathMatrixCL, reuse_expression) {
  std::string expected_kernel_src
      = "kernel void calculate(__global double* var1_global, int var1_rows, "
        "int var1_view, __global double* var2_global, int var2_rows, int "
        "var2_view, __global double* var5_global, int var5_rows, int "
        "var5_view){\n"
        "int i = get_global_id(0);\n"
        "int j = get_global_id(1);\n"
        "double var1 = 0; if (!((!contains_nonzero(var1_view, LOWER) && j < i) "
        "|| (!contains_nonzero(var1_view, UPPER) && j > i))) {var1 = "
        "var1_global[i + var1_rows * j];}\n"
        "double var2 = 0; if (!((!contains_nonzero(var2_view, LOWER) && j < i) "
        "|| (!contains_nonzero(var2_view, UPPER) && j > i))) {var2 = "
        "var2_global[i + var2_rows * j];}\n"
        "double var3 = var1+var2;\n"
        "double var4 = var3*var3;\n"
        "var5_global[i + var5_rows * j] = var4;}";
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + m2_cl;
  auto tmp2 = stan::math::elewise_multiplication(tmp, tmp);
  matrix_cl<double> res_cl;
  std::string kernel_src = tmp2.get_kernel_source_for_evaluating_into(res_cl);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  auto tmp_eig = m1.array() + m2.array();
  MatrixXd res_eig = tmp_eig * tmp_eig;

  EXPECT_MATRIX_NEAR(res_eig, res, 1e-9);
}

#endif
