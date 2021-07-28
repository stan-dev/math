#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <test/unit/math/opencl/kernel_generator/reference_kernel.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <algorithm>
#include <string>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using stan::math::matrix_cl;

TEST(KernelGenerator, addition_test) {
  std::string kernel_filename = "binary_operation_addition.cl";
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);

  auto tmp = m1_cl + m2_cl;

  matrix_cl<double> res_cl;
  std::string kernel_src = tmp.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1 + m2.cast<double>();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

#define BINARY_OPERATION_TEST(test_name, operation, res_type)          \
  TEST(KernelGenerator, test_name) {                                   \
    MatrixXd m1(3, 3);                                                 \
    m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;                           \
    MatrixXi m2(3, 3);                                                 \
    m2 << 1, 100, 1000, 0, -10, -12, 2, -8, 8;                         \
    MatrixXd m_size(3, 2);                                             \
    m_size << 1, 100, 1000, 0, -10, -12;                               \
                                                                       \
    matrix_cl<double> m1_cl(m1);                                       \
    matrix_cl<int> m2_cl(m2);                                          \
    matrix_cl<double> m_size_cl(m_size);                               \
    matrix_cl<res_type> unused_res;                                    \
                                                                       \
    EXPECT_THROW(unused_res = m1_cl operation m_size_cl,               \
                 std::invalid_argument);                               \
                                                                       \
    auto tmp = m1_cl operation m2_cl;                                  \
    matrix_cl<res_type> res_cl = tmp;                                  \
    Matrix<res_type, -1, -1> res = stan::math::from_matrix_cl(res_cl); \
                                                                       \
    Matrix<res_type, -1, -1> correct                                   \
        = m1.array() operation m2.cast<double>().array();              \
    EXPECT_MATRIX_NEAR(res, correct, 1e-9);                            \
  }

BINARY_OPERATION_TEST(subtraction_test, -, double);

TEST(KernelGenerator, elt_multiply_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);

  auto tmp = elt_multiply(m1_cl, m2_cl);
  matrix_cl<double> res_cl = tmp;
  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.array() * m2.cast<double>().array();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

TEST(KernelGenerator, elt_divide_test) {
  MatrixXd m1(3, 3);
  m1 << 1, 2.5, 3, 4, 5, 6.3, 7, -8, -9.5;
  MatrixXi m2(3, 3);
  m2 << 10, 100, 1000, 1, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<int> m2_cl(m2);
  auto tmp = elt_divide(m1_cl, m2_cl);
  matrix_cl<double> res_cl = tmp;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  MatrixXd correct = m1.array() / m2.cast<double>().array();
  EXPECT_MATRIX_NEAR(res, correct, 1e-9);
}

BINARY_OPERATION_TEST(less_than_test, <, bool);
BINARY_OPERATION_TEST(less_than_or_equal_test, <=, bool);
BINARY_OPERATION_TEST(greater_than_test, >, bool);
BINARY_OPERATION_TEST(greater_than_or_equal_test, >=, bool);
BINARY_OPERATION_TEST(equals_test, ==, bool);
BINARY_OPERATION_TEST(not_equals_test, !=, bool);

TEST(KernelGenerator, logical_or_test) {
  Matrix<bool, -1, -1> m1(3, 3);
  m1 << true, true, true, false, false, true, true, false, false;
  Matrix<bool, -1, -1> m2(3, 3);
  m2 << true, false, false, true, false, true, false, true, false;

  matrix_cl<bool> m1_cl(m1);
  matrix_cl<bool> m2_cl(m2);

  auto tmp = m1_cl || m2_cl;
  matrix_cl<bool> res_cl = tmp;
  Matrix<bool, -1, -1> res = stan::math::from_matrix_cl(res_cl);

  Matrix<bool, -1, -1> correct = m1 || m2;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, logical_and_test) {
  Matrix<bool, -1, -1> m1(3, 3);
  m1 << true, true, true, false, false, true, true, false, false;
  Matrix<bool, -1, -1> m2(3, 3);
  m2 << true, false, false, true, false, true, false, true, false;

  matrix_cl<bool> m1_cl(m1);
  matrix_cl<bool> m2_cl(m2);

  auto tmp = m1_cl && m2_cl;
  matrix_cl<bool> res_cl = tmp;
  Matrix<bool, -1, -1> res = stan::math::from_matrix_cl(res_cl);

  Matrix<bool, -1, -1> correct = m1 && m2;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(KernelGenerator, binary_operation_multiple_operations) {
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

TEST(KernelGenerator, binary_operation_multiple_operations_accepts_lvalue) {
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

TEST(KernelGenerator, multiplication_with_scalar_test) {
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

TEST(KernelGenerator, matrix_multiplication_in_expression_test) {
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

TEST(KernelGenerator, reuse_expression_simple) {
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0.1, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = stan::math::elt_divide(m1_cl, m2_cl);
  auto tmp2 = stan::math::elt_multiply(tmp, tmp);
  matrix_cl<double> res_cl;
  std::string kernel_src = tmp2.get_kernel_source_for_evaluating_into(res_cl);
  // if the expression is correctly reused, division will only occur once in the
  // kernel
  EXPECT_EQ(1, std::count(kernel_src.begin(), kernel_src.end(), '/'));
  res_cl = tmp2;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  auto tmp_eig = m1.array() / m2.array();
  MatrixXd res_eig = tmp_eig * tmp_eig;

  EXPECT_MATRIX_NEAR(res_eig, res, 1e-9);
}

// Shows subexpressions tmp and tmp2 are reused in the kernel
TEST(KernelGenerator, reuse_expression_complicated) {
  std::string kernel_filename = "binary_operation_reuse_expression.cl";
  MatrixXd m1(3, 3);
  m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2(3, 3);
  m2 << 10, 100, 1000, 0, -10, -12, 2, 4, 8;

  matrix_cl<double> m1_cl(m1);
  matrix_cl<double> m2_cl(m2);
  auto tmp = m1_cl + m2_cl;
  auto tmp2 = stan::math::elt_divide(stan::math::elt_multiply(tmp, tmp), m1_cl);
  auto tmp3 = stan::math::elt_multiply(stan::math::elt_divide(tmp, tmp2), tmp2);
  matrix_cl<double> res_cl;
  std::string kernel_src = tmp3.get_kernel_source_for_evaluating_into(res_cl);
  stan::test::store_reference_kernel_if_needed(kernel_filename, kernel_src);
  std::string expected_kernel_src
      = stan::test::load_reference_kernel(kernel_filename);
  EXPECT_EQ(expected_kernel_src, kernel_src);

  res_cl = tmp3;

  MatrixXd res = stan::math::from_matrix_cl(res_cl);

  auto tmp_eig = m1.array() + m2.array();
  MatrixXd tmp2_eig = (tmp_eig * tmp_eig).array() / m1.array();
  MatrixXd tmp3_eig
      = (tmp_eig.array() / tmp2_eig.array()).array() * tmp2_eig.array();

  EXPECT_MATRIX_NEAR(tmp3_eig, res, 1e-9);
}

#endif
