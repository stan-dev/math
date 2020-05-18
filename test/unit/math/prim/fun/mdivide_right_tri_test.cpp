#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#include <boost/random/mersenne_twister.hpp>
#endif

TEST(MathMatrixPrim, mdivide_right_tri_val) {
  using stan::math::mdivide_right_tri;
  stan::math::matrix_d I = Eigen::MatrixXd::Identity(2, 2);

  stan::math::matrix_d Ad(2, 2);
  Ad << 2.0, 0.0, 5.0, 7.0;
  expect_matrix_eq(I, mdivide_right_tri<Eigen::Lower>(Ad, Ad));

  Ad << 2.0, 3.0, 0.0, 7.0;
  expect_matrix_eq(I, mdivide_right_tri<Eigen::Upper>(Ad, Ad));
}

TEST(MathMatrixPrim, mdivide_right_tri_size_zero) {
  using stan::math::mdivide_right_tri;
  stan::math::matrix_d Ad(0, 0);
  stan::math::matrix_d b0(2, 0);
  stan::math::matrix_d I;

  I = mdivide_right_tri<Eigen::Lower>(Ad, Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_right_tri<Eigen::Upper>(Ad, Ad);
  EXPECT_EQ(0, I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_right_tri<Eigen::Lower>(b0, Ad);
  EXPECT_EQ(b0.rows(), I.rows());
  EXPECT_EQ(0, I.cols());

  I = mdivide_right_tri<Eigen::Upper>(b0, Ad);
  EXPECT_EQ(b0.rows(), I.rows());
  EXPECT_EQ(0, I.cols());
}

#ifdef STAN_OPENCL

void mdivide_right_tri_cl_test(int size) {
  boost::random::mt19937 rng;
  stan::math::matrix_d m1(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      m1(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    m1(i, i) = 20.0;
    for (int j = i + 1; j < size; j++) {
      m1(i, j) = 0.0;
    }
  }

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = size * 2;

  stan::math::matrix_d m1_cpu
      = stan::math::mdivide_right_tri<Eigen::Lower>(m1, m1);

  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer = 0;

  stan::math::matrix_d m1_cl
      = stan::math::mdivide_right_tri<Eigen::Lower>(m1, m1);

  EXPECT_MATRIX_NEAR(m1_cpu, m1_cl, 1E-8);
}
TEST(MathMatrixCL, mdivide_right_tri_cl_small) { mdivide_right_tri_cl_test(3); }
TEST(MathMatrixCL, mdivide_right_tri_cl_mid) { mdivide_right_tri_cl_test(100); }
TEST(MathMatrixCL, mdivide_right_tri_cl_big) { mdivide_right_tri_cl_test(500); }
#endif
