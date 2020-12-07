#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathMatrixCL, inverse_cl_exception) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  stan::math::matrix_cl<double> m2(m1, stan::math::matrix_cl_view::Lower);
  stan::math::matrix_cl<double> m3(2, 3);
  using stan::math::tri_inverse;
  EXPECT_THROW(m3 = tri_inverse(m2), std::invalid_argument);
  m2.view(stan::math::matrix_cl_view::Upper);
  EXPECT_THROW(m3 = tri_inverse(m2), std::invalid_argument);

  stan::math::matrix_d m4(3, 3);
  m4 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::matrix_cl<double> m5(m1, stan::math::matrix_cl_view::Entire);
  EXPECT_THROW(m3 = tri_inverse(m5), std::invalid_argument);
  m5.view(stan::math::matrix_cl_view::Diagonal);
  EXPECT_THROW(m3 = tri_inverse(m5), std::invalid_argument);
}

void lower_inverse_test(int size) {
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

  m1_cpu = stan::math::mdivide_left_tri<Eigen::Lower>(m1);

  stan::math::matrix_cl<double> m2(m1, stan::math::matrix_cl_view::Lower);
  auto m3 = stan::math::tri_inverse(m2);
  m1_cl = stan::math::from_matrix_cl(m3);
  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j <= i; j++) {
      double abs_err = std::fabs(m1_cpu(i, j) - m1_cl(i, j));
      double a = std::max(abs_err / m1_cpu(i, j), abs_err / m1_cl(i, j));
      max_error = std::max(max_error, a);
    }
  }
  EXPECT_LT(max_error, 1e-8);
}

void upper_inverse_test(int size) {
  boost::random::mt19937 rng;
  auto m1 = stan::math::matrix_d(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      m1(i, j) = 0.0;
    }
    m1(i, i) = 10000.0;
    for (int j = i + 1; j < size; j++) {
      m1(i, j) = stan::math::uniform_rng(-5, 5, rng);
    }
  }
  stan::math::matrix_d m1_cpu(size, size);
  stan::math::matrix_d m1_cl(size, size);

  m1_cpu = stan::math::mdivide_left_tri<Eigen::Upper>(m1);

  stan::math::matrix_cl<double> m2(m1, stan::math::matrix_cl_view::Upper);
  auto m3 = stan::math::tri_inverse(m2);
  m1_cl = stan::math::from_matrix_cl(m3);
  double max_error = 0;
  for (int i = 0; i < size; i++) {
    for (int j = i; j < size; j++) {
      double abs_err = std::fabs(m1_cpu(i, j) - m1_cl(i, j));
      double a = std::max(abs_err / m1_cpu(i, j), abs_err / m1_cl(i, j));
      max_error = std::max(max_error, a);
    }
  }
  EXPECT_LT(max_error, 1e-8);
}

TEST(MathMatrixCL, inverse_cl_small) {
  lower_inverse_test(3);
  upper_inverse_test(3);
}

TEST(MathMatrixCL, inverse_cl_1under_block_size) {
  lower_inverse_test(31);
  upper_inverse_test(31);
}

TEST(MathMatrixCL, inverse_cl_1over_block_size) {
  lower_inverse_test(33);
  upper_inverse_test(33);
}

TEST(MathMatrixCL, inverse_cl_big_power_of_2) {
  lower_inverse_test(512);
  upper_inverse_test(512);
}

TEST(MathMatrixCL, inverse_cl_big_non_power_of_2) {
  lower_inverse_test(700);
  upper_inverse_test(700);
}

TEST(MathMatrixCL, inverse_cl_very_big_non_power_of_2) {
  lower_inverse_test(1500);
  upper_inverse_test(1500);
}
#endif
