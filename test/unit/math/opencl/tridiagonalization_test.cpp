#ifdef STAN_OPENCL
#include <stan/math/opencl/tridiagonalization.hpp>
#include <stan/math/opencl/prim/identity_matrix.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixGPU, tridiagonalization_trivial) {
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(3, 3);
  stan::math::matrix_cl<double> id_cl(id);
  stan::math::matrix_cl<double> packed_cl;
  stan::math::internal::block_householder_tridiag_cl(id_cl, packed_cl);
  EXPECT_NEAR_REL(stan::math::from_matrix_cl(packed_cl),
                  stan::math::from_matrix_cl(id_cl));
  stan::math::matrix_cl<double> q_cl = id_cl;
  stan::math::internal::block_apply_packed_Q_cl(packed_cl, q_cl);
  EXPECT_NEAR_REL(stan::math::from_matrix_cl(q_cl),
                  stan::math::from_matrix_cl(id_cl));
}

TEST(MathMatrixGPU, tridiagonalization_small) {
  int size = 7;
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size, size);
  input += input.transpose().eval();
  stan::math::matrix_cl<double> input_cl(input);
  stan::math::matrix_cl<double> packed_cl;

  stan::math::internal::block_householder_tridiag_cl(input_cl, packed_cl);
  stan::math::matrix_cl<double> q_cl
      = stan::math::identity_matrix<stan::math::matrix_cl<double>>(size);
  Eigen::MatrixXd q2_in = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd q2 = q2_in;
  stan::math::matrix_cl<double> q2_cl(q2);
  stan::math::internal::block_apply_packed_Q_cl(packed_cl, q_cl);
  stan::math::internal::block_apply_packed_Q_cl(packed_cl, q2_cl);

  Eigen::MatrixXd packed = stan::math::from_matrix_cl(packed_cl);
  Eigen::MatrixXd q = stan::math::from_matrix_cl(q_cl);
  Eigen::MatrixXd q2_res = stan::math::from_matrix_cl(q2_cl);
  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = packed.diagonal();
  t.diagonal(1) = packed.diagonal(1);
  t.diagonal(-1) = packed.diagonal(1);
  EXPECT_TRUE(
      (q * q.transpose()).isApprox(Eigen::MatrixXd::Identity(size, size)));
  EXPECT_TRUE((q * t * q.transpose()).isApprox(input));
  EXPECT_TRUE((q * q2_in).isApprox(q2_res));
}

TEST(MathMatrixGPU, tridiagonalization_large) {
  int size = 345;
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size, size);
  input += input.transpose().eval();
  stan::math::matrix_cl<double> input_cl(input);
  stan::math::matrix_cl<double> packed_cl;

  stan::math::internal::block_householder_tridiag_cl(input_cl, packed_cl);
  stan::math::matrix_cl<double> q_cl
      = stan::math::identity_matrix<stan::math::matrix_cl<double>>(size);
  Eigen::MatrixXd q2_in = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd q2 = q2_in;
  stan::math::matrix_cl<double> q2_cl(q2);
  stan::math::matrix_cl<double> q2_2_cl(q2_cl);
  stan::math::matrix_cl<double> q1_2_cl(q_cl);
  stan::math::internal::block_apply_packed_Q_cl(packed_cl, q_cl);
  stan::math::internal::block_apply_packed_Q_cl(packed_cl, q2_cl);

  Eigen::MatrixXd packed = stan::math::from_matrix_cl(packed_cl);
  Eigen::MatrixXd q = stan::math::from_matrix_cl(q_cl);
  Eigen::MatrixXd q2_res = stan::math::from_matrix_cl(q2_cl);
  Eigen::MatrixXd t = Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal() = packed.diagonal();
  t.diagonal(1) = packed.diagonal(1);
  t.diagonal(-1) = packed.diagonal(1);
  EXPECT_TRUE(
      (q * q.transpose()).isApprox(Eigen::MatrixXd::Identity(size, size)));
  EXPECT_TRUE((q * t * q.transpose()).isApprox(input));
  EXPECT_TRUE((q * q2_in).isApprox(q2_res));
}
#endif
