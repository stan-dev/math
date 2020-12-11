#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

auto matrix_multiply_functor
    = [](const auto& a, const auto& b) { return stan::math::multiply(a, b); };

TEST(OpenCLMatrixMultiply, prim_rev_values_small) {
  int N = 2;
  int M = 3;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(M, K);
  b << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_N_0) {
  int N = 0;
  int M = 3;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(M, K);
  b << 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1;
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_M_0) {
  int N = 2;
  int M = 0;
  int K = 4;

  Eigen::MatrixXd a(N, M);
  Eigen::MatrixXd b(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_K_0) {
  int N = 2;
  int M = 3;
  int K = 0;

  Eigen::MatrixXd a(N, M);
  a << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_values_large) {
  int N = 71;
  int M = 83;
  int K = 97;

  Eigen::MatrixXd a = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(M, K);
  stan::math::test::compare_cpu_opencl_prim_rev(matrix_multiply_functor, a, b);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_zero) {
  int N = 0;
  int M = 0;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<Eigen::Matrix<double, -1, -1>> A1 = A;
  var_value<double> a1 = 2.0;
  var_value<Eigen::Matrix<double, -1, -1>> B1 = stan::math::multiply(a1, A1);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  var_value<matrix_cl<double>> A2_cl = stan::math::to_matrix_cl(A2);
  var_value<double> a2 = 2.0;
  var_value<matrix_cl<double>> B2_cl = stan::math::multiply(a2, A2_cl);
  stan::math::var C2 = stan::math::sum(B2_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B2_cl.val()), 1E-8);
  EXPECT_MATRIX_NEAR(A1.adj(), A2.adj(), 1E-8);
  EXPECT_NEAR(a1.adj(), a2.adj(), 1E-8);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_mat) {
  int N = 71;
  int M = 83;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<Eigen::Matrix<double, -1, -1>> A1 = A;
  var_value<double> a1 = 2.0;
  var_value<Eigen::Matrix<double, -1, -1>> B1 = stan::math::multiply(a1, A1);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  var_value<matrix_cl<double>> A2_cl = stan::math::to_matrix_cl(A2);
  var_value<double> a2 = 2.0;
  var_value<matrix_cl<double>> B2_cl = stan::math::multiply(a2, A2_cl);
  stan::math::var C2 = stan::math::sum(B2_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B2_cl.val()), 1E-8);
  EXPECT_MATRIX_NEAR(A1.adj(), A2.adj(), 1E-8);
  EXPECT_NEAR(a1.adj(), a2.adj(), 1E-8);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_d_mat_v) {
  int N = 71;
  int M = 83;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<Eigen::Matrix<double, -1, -1>> A1 = A;
  double a1 = 2.0;
  var_value<Eigen::Matrix<double, -1, -1>> B1 = stan::math::multiply(a1, A1);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  var_value<matrix_cl<double>> A2_cl = stan::math::to_matrix_cl(A2);
  double a2 = 2.0;
  var_value<matrix_cl<double>> B2_cl = stan::math::multiply(a2, A2_cl);
  stan::math::var C2 = stan::math::sum(B2_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B2_cl.val()), 1E-8);
  EXPECT_MATRIX_NEAR(A1.adj(), A2.adj(), 1E-8);
}

TEST(OpenCLMatrixMultiply, prim_rev_scalar_v_mat_d) {
  int N = 71;
  int M = 83;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<double> a1 = 2.0;
  Eigen::Matrix<stan::math::var, -1, -1> B1 = stan::math::multiply(a1, A);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  matrix_cl<double> A_cl = stan::math::to_matrix_cl(A);
  var_value<double> a2 = 2.0;
  var_value<matrix_cl<double>> B_cl = stan::math::multiply(a2, A_cl);
  stan::math::var C2 = stan::math::sum(B_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B_cl.val()), 1E-8);
  EXPECT_NEAR(a1.adj(), a2.adj(), 1E-8);
}

TEST(OpenCLMatrixMultiply, prim_rev_mat_scalar_zero) {
  int N = 0;
  int M = 0;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<Eigen::Matrix<double, -1, -1>> A1 = A;
  var_value<double> a1 = 4.0;
  var_value<Eigen::Matrix<double, -1, -1>> B1 = stan::math::multiply(A1, a1);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  var_value<matrix_cl<double>> A2_cl = stan::math::to_matrix_cl(A2);
  var_value<double> a2 = 4.0;
  var_value<matrix_cl<double>> B2_cl = stan::math::multiply(A2_cl, a2);
  stan::math::var C2 = stan::math::sum(B2_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B2_cl.val()), 1E-8);
  EXPECT_MATRIX_NEAR(A1.adj(), A2.adj(), 1E-8);
  EXPECT_NEAR(a1.adj(), a2.adj(), 1E-8);
}

TEST(OpenCLMatrixMultiply, prim_rev_mat_scalar) {
  int N = 71;
  int M = 83;
  using stan::math::var_value;
  using stan::math::matrix_cl;
  Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);

  var_value<Eigen::Matrix<double, -1, -1>> A1 = A;
  var_value<double> a1 = 4.0;
  var_value<Eigen::Matrix<double, -1, -1>> B1 = stan::math::multiply(A1, a1);
  stan::math::var C1 = stan::math::sum(B1);

  var_value<Eigen::Matrix<double, -1, -1>> A2 = A;
  var_value<matrix_cl<double>> A2_cl = stan::math::to_matrix_cl(A2);
  var_value<double> a2 = 4.0;
  var_value<matrix_cl<double>> B2_cl = stan::math::multiply(A2_cl, a2);
  stan::math::var C2 = stan::math::sum(B2_cl);
  (C1 + C2).grad();

  EXPECT_MATRIX_NEAR(B1.val(), stan::math::from_matrix_cl(B2_cl.val()), 1E-8);
  EXPECT_MATRIX_NEAR(A1.adj(), A2.adj(), 1E-8);
  EXPECT_NEAR(a1.adj(), a2.adj(), 1E-8);
}


#endif
