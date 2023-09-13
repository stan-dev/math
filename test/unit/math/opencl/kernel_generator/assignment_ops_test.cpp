#ifdef STAN_OPENCL
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <string>

TEST(KernelGenerator, plus_equals) {
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(10, 10);
  matrix_cl<double> A_cl = to_matrix_cl(A);
  matrix_cl<double> B_cl = to_matrix_cl(B);
  matrix_cl<double> C_cl = to_matrix_cl(C);
  C += A + B;
  results(C_cl) += expressions(A_cl + B_cl);
  Eigen::MatrixXd C_cl_host = from_matrix_cl(C_cl);
  EXPECT_MATRIX_EQ(C_cl_host, C)
}

TEST(KernelGenerator, minus_equals) {
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(10, 10);
  matrix_cl<double> A_cl = to_matrix_cl(A);
  matrix_cl<double> B_cl = to_matrix_cl(B);
  matrix_cl<double> C_cl = to_matrix_cl(C);
  C -= A + B;
  results(C_cl) -= expressions(A_cl + B_cl);
  Eigen::MatrixXd C_cl_host = from_matrix_cl(C_cl);
  EXPECT_MATRIX_EQ(C_cl_host, C)
}

TEST(KernelGenerator, divide_equals) {
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(10, 10);
  matrix_cl<double> A_cl = to_matrix_cl(A);
  matrix_cl<double> B_cl = to_matrix_cl(B);
  matrix_cl<double> C_cl = to_matrix_cl(C);
  C.array() /= A.array() + B.array();
  results(C_cl) /= expressions(A_cl + B_cl);
  Eigen::MatrixXd C_cl_host = from_matrix_cl(C_cl);
  EXPECT_MATRIX_EQ(C_cl_host, C)
}

TEST(KernelGenerator, times_equals) {
  using stan::math::from_matrix_cl;
  using stan::math::matrix_cl;
  using stan::math::to_matrix_cl;
  using stan::math::var;
  using stan::math::var_value;

  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd C = Eigen::MatrixXd::Random(10, 10);
  matrix_cl<double> A_cl = to_matrix_cl(A);
  matrix_cl<double> B_cl = to_matrix_cl(B);
  matrix_cl<double> C_cl = to_matrix_cl(C);
  C.array() *= A.array() + B.array();
  results(C_cl) *= expressions(A_cl + B_cl);
  Eigen::MatrixXd C_cl_host = from_matrix_cl(C_cl);
  EXPECT_MATRIX_EQ(C_cl_host, C)
}

#endif
