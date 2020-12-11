#ifdef STAN_OPENCL

#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(MathGpu, make_kernel) {
  stan::math::matrix_d m(3, 3);

  stan::math::matrix_cl<double> m_cl(m.cols(), m.rows());
  stan::math::opencl_kernels::fill_strict_tri(
      cl::NDRange(m_cl.rows(), m_cl.cols()), m_cl, 0, m_cl.rows(), m_cl.cols(),
      stan::math::matrix_cl_view::Entire);
  m = stan::math::from_matrix_cl(m_cl);
}

TEST(MathGpu, write_after_write) {
  using stan::math::matrix_cl;

  for (int j = 0; j < 1000; j++) {
    matrix_cl<double> m_cl(3, 3);

    for (int i = 0; i < 4; i++) {
      stan::math::opencl_kernels::fill_strict_tri(
          cl::NDRange(3, 3), m_cl, i, 3, 3, stan::math::matrix_cl_view::Entire);
    }

    stan::math::matrix_d res = stan::math::from_matrix_cl(m_cl);
    stan::math::matrix_d correct = stan::math::matrix_d::Constant(3, 3, 3);

    EXPECT_MATRIX_NEAR(res, correct, 1e-13);
  }
}

TEST(MathGpu, read_after_write_and_write_after_read) {
  using stan::math::matrix_cl;

  matrix_cl<double> m_cl(3, 3);
  matrix_cl<double> ms_cl[1000];
  for (int j = 0; j < 1000; j++) {
    stan::math::opencl_kernels::fill_strict_tri(
        cl::NDRange(3, 3), m_cl, j, 3, 3, stan::math::matrix_cl_view::Entire);

    ms_cl[j] = m_cl;
  }
  for (int j = 0; j < 1000; j++) {
    stan::math::matrix_d res = stan::math::from_matrix_cl(ms_cl[j]);
    stan::math::matrix_d correct = stan::math::matrix_d::Constant(3, 3, j);

    EXPECT_MATRIX_NEAR(res, correct, 1e-13);
  }
}

#endif
