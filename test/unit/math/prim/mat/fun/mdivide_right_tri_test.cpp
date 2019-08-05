#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#include <boost/random/mersenne_twister.hpp>
#endif

#define EXPECT_MATRIX_NEAR(A, B, DELTA) \
  for (int i = 0; i < A.size(); i++)    \
    EXPECT_NEAR(A(i), B(i), DELTA);

TEST(MathMatrix, mdivide_right_tri_val) {
  using stan::math::mdivide_right_tri;
  stan::math::matrix_d Ad(2, 2);
  stan::math::matrix_d I;

  Ad << 2.0, 0.0, 5.0, 7.0;

  I = mdivide_right_tri<Eigen::Lower>(Ad, Ad);
  EXPECT_NEAR(1.0, I(0, 0), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1), 1.0e-12);

  Ad << 2.0, 3.0, 0.0, 7.0;

  I = mdivide_right_tri<Eigen::Upper>(Ad, Ad);
  EXPECT_NEAR(1.0, I(0, 0), 1.0E-12);
  EXPECT_NEAR(0.0, I(0, 1), 1.0E-12);
  EXPECT_NEAR(0.0, I(1, 0), 1.0E-12);
  EXPECT_NEAR(1.0, I(1, 1), 1.0e-12);
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
