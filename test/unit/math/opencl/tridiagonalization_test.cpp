#ifdef STAN_OPENCL
#include <stan/math/opencl/tridiagonalization.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixGPU, tridiagonalization_trivial) {
  Eigen::MatrixXd id = Eigen::MatrixXd::Identity(3,3);
  Eigen::MatrixXd packed;
  stan::math::internal::block_householder_tridiag_cl(id,packed);
  EXPECT_TRUE(packed.isApprox(id));
  Eigen::MatrixXd q = id;
  stan::math::internal::block_apply_packed_Q_cl(packed,q);
  EXPECT_TRUE(q.isApprox(id));
}

TEST(MathMatrixGPU, tridiagonalization_small) {
  int size=7;
  srand(0); //ensure test repeatability
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size,size);
  input += input.transpose().eval();
  Eigen::MatrixXd packed;

  stan::math::internal::block_householder_tridiag_cl(input,packed);
  Eigen::MatrixXd q = Eigen::MatrixXd::Identity(size,size);
  stan::math::internal::block_apply_packed_Q_cl(packed,q);

  Eigen::MatrixXd t=Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal()=packed.diagonal();
  t.diagonal(1)=packed.diagonal(1);
  t.diagonal(-1)=packed.diagonal(1);
  EXPECT_TRUE((q*q.transpose()).isApprox(Eigen::MatrixXd::Identity(size,size)));
  EXPECT_TRUE((q*t*q.transpose()).isApprox(input));
}

TEST(MathMatrixGPU, tridiagonalization_large) {
  int size=345;
  srand(0); //ensure test repeatability
  Eigen::MatrixXd input = Eigen::MatrixXd::Random(size,size);
  input += input.transpose().eval();
  Eigen::MatrixXd packed;

  stan::math::internal::block_householder_tridiag_cl(input,packed);
  Eigen::MatrixXd q = Eigen::MatrixXd::Identity(size,size);
  stan::math::internal::block_apply_packed_Q_cl(packed,q);

  Eigen::MatrixXd t=Eigen::MatrixXd::Constant(size, size, 0);
  t.diagonal()=packed.diagonal();
  t.diagonal(1)=packed.diagonal(1);
  t.diagonal(-1)=packed.diagonal(1);
  EXPECT_TRUE((q*q.transpose()).isApprox(Eigen::MatrixXd::Identity(size,size),1e-9));
  EXPECT_TRUE((q*t*q.transpose()).isApprox(input,1e-9));
}
#endif