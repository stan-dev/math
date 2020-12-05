#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixCL, mdivide_left_tri_low_cl_exception) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> m2(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_cl<double> m3(2, 3);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m2), std::invalid_argument);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m2, m2), std::invalid_argument);
  m2.view(stan::math::matrix_cl_view::Upper);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m2), std::invalid_argument);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m2, m2), std::invalid_argument);

  stan::math::matrix_d m4(3, 3);
  m4 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl<double> m5(m1, stan::math::matrix_cl_view::Entire);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m5), std::invalid_argument);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m5, m5), std::invalid_argument);
  m5.view(stan::math::matrix_cl_view::Diagonal);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m5), std::invalid_argument);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(m5, m5), std::invalid_argument);
}

void mdivide_left_tri_low_Ab_test(int size) {
  boost::random::mt19937 rng;
  auto m1 = stan::math::matrix_d(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      m1(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    m1(i, i) = 10000.0;
    for (int j = i + 1; j < size; j++) {
      m1(i, j) = 0.0;
    }
    m1(i, 0) = stan::math::uniform_rng(-5, 5, rng);
  }
  // force the CPU version for comparison
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = size * 2;
  stan::math::matrix_d m1_cpu = stan::math::mdivide_left_tri_low(m1, m1);

  stan::math::matrix_cl<double> m0_cl(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_d m1_cl = stan::math::from_matrix_cl(
      stan::math::mdivide_left_tri_low(m0_cl, m0_cl));
  double max_error = 0;
  for (int i = 0; i < m1_cpu.size(); i++) {
    stan::test::expect_near_rel("mdivide_left_tri_low (OpenCL)", m1_cpu(i),
                                m1_cl(i));
  }
}

void mdivide_left_tri_low_A_test(int size) {
  boost::random::mt19937 rng;
  auto m1 = stan::math::matrix_d(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      m1(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
    m1(i, i) = 10000.0;
    for (int j = i + 1; j < size; j++) {
      m1(i, j) = 0.0;
    }
  }
  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);
  // force the CPU version for comparison
  stan::math::opencl_context.tuning_opts().tri_inverse_size_worth_transfer
      = size * 2;
  m1_cpu = stan::math::mdivide_left_tri_low(m1);

  stan::math::matrix_cl<double> m1_A_cl(m1, stan::math::matrix_cl_view::Lower);
  m1_cl = stan::math::from_matrix_cl(stan::math::mdivide_left_tri_low(m1_A_cl));
  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j <= i; j++) {
      stan::test::expect_near_rel("mdivide_left_tri_low (OpenCL)", m1_cpu(i, j),
                                  m1_cl(i, j));
    }
  }
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_small) {
  mdivide_left_tri_low_Ab_test(3);
  mdivide_left_tri_low_A_test(3);
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_1under_block_size) {
  mdivide_left_tri_low_Ab_test(31);
  mdivide_left_tri_low_A_test(31);
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_1over_block_size) {
  mdivide_left_tri_low_Ab_test(33);
  mdivide_left_tri_low_A_test(33);
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_big_power_of_2) {
  mdivide_left_tri_low_Ab_test(512);
  mdivide_left_tri_low_A_test(512);
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_big_non_power_of_2) {
  mdivide_left_tri_low_Ab_test(700);
  mdivide_left_tri_low_A_test(700);
}

TEST(MathMatrixCL, mdivide_left_tri_low_test_cl_very_big_non_power_of_2) {
  mdivide_left_tri_low_Ab_test(1500);
  mdivide_left_tri_low_A_test(1500);
}
#endif
