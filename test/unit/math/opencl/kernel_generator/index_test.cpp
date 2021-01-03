#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(KernelGenerator, row_index_errors_test) {
  EXPECT_THROW(stan::math::row_index().eval(), std::invalid_argument);
  EXPECT_NO_THROW(stan::math::row_index(2, 3).eval());
}

TEST(KernelGenerator, col_index_errors_test) {
  EXPECT_THROW(stan::math::col_index().eval(), std::invalid_argument);
  EXPECT_NO_THROW(stan::math::col_index(2, 3).eval());
}

TEST(KernelGenerator, row_index_test) {
  stan::math::matrix_cl<int> m_cl = stan::math::row_index(3, 5);

  Eigen::MatrixXi correct(3, 5);
  correct << 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2;
  Eigen::MatrixXi res1 = stan::math::from_matrix_cl(m_cl);
  EXPECT_MATRIX_EQ(res1, correct);

  m_cl = m_cl + stan::math::row_index();

  Eigen::MatrixXi res2 = stan::math::from_matrix_cl(m_cl);
  EXPECT_MATRIX_EQ(res2, correct * 2);
}

TEST(KernelGenerator, col_index_test) {
  stan::math::matrix_cl<int> m_cl = stan::math::col_index(3, 5);

  Eigen::MatrixXi correct(3, 5);
  correct << 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4;
  Eigen::MatrixXi res1 = stan::math::from_matrix_cl(m_cl);
  EXPECT_MATRIX_EQ(res1, correct);

  m_cl = m_cl + stan::math::col_index();

  Eigen::MatrixXi res2 = stan::math::from_matrix_cl(m_cl);
  EXPECT_MATRIX_EQ(res2, correct * 2);
}

#endif
